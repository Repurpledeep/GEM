package Identify_Geno_functions;
use strict;
use warnings;

use File::Basename;

use identify::Fasta_deal;
use identify::Run_jobs;
use identify::Configure;
use Exporter;

my %Bin;
my ($p_num,$runtype,$queue,$retain);

sub Geno_Workflow
{
    my $opt = shift;
    Configure::Save_opt($opt);
    ($p_num,$runtype,$queue,$retain) = Configure::Load_opt('Process','Runtype','Queue','Retain');
    %Bin = Configure::Bin_Adr();
    &Make_dir($$opt{'Outdir'});
    my $tblastn_out = &Input_VS_genome($opt);
    my $tbn_fa = &Extract_fa_from_genome($tblastn_out,$opt);
    my $blastx_out = &Candidate_VS_validating_DB($tbn_fa,$opt);
    my %fa_pep_pair = &Validated_fa_result($blastx_out,$tbn_fa,$opt);
    my $refine_fa = &Genewise_flow(\%fa_pep_pair,$opt);
}

sub Make_dir
{
    my $out = shift;
    mkdir "$out" unless(-e $out);
    mkdir "$out/tblastn" unless(-e "$out/tblastn");
    mkdir "$out/blastx" unless(-e "$out/blastx");
    mkdir "$out/genewise" unless (-e "$out/genewise");
    mkdir "$out/genewise/data" unless(-e "$out/genewise/data");

}

sub Input_VS_genome
{
    my $opt = shift;
    my ($pep,$ref,$tbn_eval,$out_dir) = ($$opt{'Input'},$$opt{'Ref'},$$opt{'T_Eval'},"$$opt{'Outdir'}/tblastn");

    my $genome_name = basename $ref;
    die "Please make sure the name of Ref is end with .fa\n" unless($genome_name =~s/\.fa$//);
    
    my @sub_db = Fasta_deal::Split_fasta_average_size($p_num,$ref,$out_dir); 
    my @subfile;
    push @subfile,"$pep";
    my %Param;
    %{$Param{'tblastn'}} = (
            '-evalue' => $tbn_eval, 
#            '-max_target_seqs' => 1,
            '-outfmt' => 6,
            );
#   my $conf = join " ",%{$Param->{'blastall'}};
    my $work_sh = "$out_dir/$genome_name\_tbn.sh";
    my @output = Run_jobs::Create_blastall(\@subfile,\@sub_db,\%Param,$work_sh,$out_dir);
    Run_jobs::Run_jobs($work_sh,$runtype,$p_num,$queue);
    my $merge_file = Run_jobs::Merge_file(\@output,"$genome_name\_tblastn",$out_dir);
    return "$merge_file";
}

sub Extract_fa_from_genome
{
    my ($tbn,$opt) = @_;
    my $ref = $$opt{'Ref'};
    my ($gff,$fasta) = ($$opt{'Gff'},$$opt{'Input'}) if($$opt{'Gff'});
###### Load and filter
    my @tbn_out = &Load_tbn_result($tbn,$gff,$fasta);
#    print "totally ",int(@tbn_out)," tbn hints\n";
#    exit 1;
###### Merge result
    my %regions = &Merge_regions_result(\@tbn_out,2,9,10);
    my @chr =keys %regions;
    my $count = 0;
    foreach(@chr)
    {
        $count += int(@{$regions{$_}});
    }
    print "Merge confirm:$count\n";
    Run_jobs::Print_time();
#    exit 1;
    my $output = "$tbn.fa";
###### Extract
    &Check_with_gff(\%regions,$fasta,$gff) if($gff);   ###$regions{$chr}->@regions
    Fasta_deal::Extract_by_regions(\%regions,$ref,$output);
    return $output;
#    exit 1;
}

sub Load_tbn_result             
{
    my ($tbn,$gff,$fasta) = @_;
    my ($aa_len,$identity,$max_gap) = Configure::Load_opt('AA_len','Identity','Max_gap');
    print "Load_gff\n";
    Run_jobs::Print_time();
    my %gff = Fasta_deal::Load_gff_byID($gff,$fasta) if($gff);  # %gff Structure : $gff{chr} = \@regions
    my @output;   #### whole orginal line
    my $total;
    print "Load_tbn\n";
    Run_jobs::Print_time();
    open TBN,"<$tbn" or die "$!";
    while(<TBN>)
    {
        chomp;
        my @word = split /\t/,$_;
        my $identity_len = $word[2] * $word[3] / 100;    # amino acid length; identity is percent, so identity/100.
        $total++;
        if($aa_len && $identity_len < $aa_len)
        {
#            print "AA\n";
            next;
        }
        elsif($identity && $word[2] < $identity)
        {
#            print "Iden\n";
           next;
        }
        elsif($max_gap && $word[5] > $max_gap)
        {
#            print "gap\n";
            next;
        }
        else
        {
            if($gff)
            {
                my $strand = ($word[9] > $word[8]) ? '+' : '-';
                my @tbn_region = ($strand eq '+') ? ($word[1],$word[8],$word[9],$strand) : ($word[1],$word[9],$word[8],$strand);
                my $judge = &Check_with_gff_Ver2(\@{$gff{$word[1]}},\@tbn_region);
                if($judge == 1)
                {
#                    print "Function check\n";
                    next;
                }
                else
                {
                    push @output,$_;
                }
            }
            else
            {
                push @output,$_;
            }
        }
    }
    close TBN;
    print "Finished Loading:$total\n";
    Run_jobs::Print_time();
    return @output;
}

sub Check_with_gff_Ver2
{
    my ($gff,$tbn) = @_;
    my $max_overlap = Configure::Load_opt('Overlap');
#    print "$max_overlap\n";
    for(my $i = 0;$i < @{$gff}; $i++)
    {
        my @gff_region = split /\t/,$$gff[$i];
        my $gene_len = abs($gff_region[2] - $gff_region[1]) + 1;
        my $merge = &Compare_region($tbn,\@gff_region,'merge');
        if($merge == 1)
        {
            my $new_len = abs($gff_region[2] - $gff_region[1]) + 1;
            my $tbn_len = abs($$tbn[2] - $$tbn[1]) + 1;
            my $overlap_len = $tbn_len - ($new_len - $gene_len);
            my $overlap_ratio = $overlap_len / $tbn_len;
#            print "Function_Check:$$gff[$i]\t$$tbn[1]\t$$tbn[2]\t$overlap_len\t$tbn_len\t$overlap_ratio\n";
            if($overlap_ratio > $max_overlap)
            {
#                print "Check_ver2:$$gff[$i]\t$$tbn[1]\t$$tbn[2]\t$overlap_len\t$tbn_lent\t$overlap_ratio\t$max_overlap\n";
                return 1;
            }
            else
            {
                next;
            }
        }
    }
    return 0;
}

sub Candidate_VS_validating_DB
{
    my ($tbn_fa,$opt) = @_;
    my ($p_num,$eval) = ($$opt{'Process'},$$opt{'P_Eval'});
    my $out_dir = "$$opt{'Outdir'}/blastx";
    my $basename = basename $tbn_fa;
    
    my @subfiles = Fasta_deal::Split_fasta_average_size($p_num,$tbn_fa,$out_dir);
    my @sub_db; push @sub_db,$$opt{'DB'};
    my %Param;
    %{$Param{'blastx'}} = (
            '-evalue' => $eval,
            '-max_target_seqs' => 1,
            '-strand' => 'plus',
            '-outfmt' => 6,
            );
#    my $conf= join (" ",%{$Param->{'blastx'}});
    my $shell_file = "$out_dir/$basename\_blastx.sh";

    my @output = Run_jobs::Create_blastall(\@subfiles,\@sub_db,\%Param,$shell_file,$out_dir);
    Run_jobs::Run_jobs($shell_file,$runtype,$p_num,$queue);
    my $merge_file = Run_jobs::Merge_file(\@output,"$basename\_blastx",$out_dir);
    return $merge_file;
}

sub Validated_fa_result
{
    my ($blastx,$tbn_fa,$opt) = @_;
    my ($outdir,$DB) = ("$$opt{'Outdir'}/genewise/data/",$$opt{'DB'});
    my %tbn = Fasta_deal::Load_fasta_by_id($tbn_fa);
    my %db = Fasta_deal::Load_fasta_by_id($DB);
    my %pairs;
    my %uniq;
    my $count = 1;
    open BLAX,"<$blastx" or die "$!";
    while(<BLAX>)
    {
        chomp;
        my @word = split /\t/;
        next unless($word[1] =~/NBS-LRR/);
        if(!$tbn{$word[0]} or !$db{$word[1]})
        {
            die "Can't find ID : $word[0] or $word[1]\n";
        }
        else
        {
            unless($uniq{$word[0]})
            {
                open DNA,">$outdir/$count.fa" or die "$!";
                print DNA ">$word[0]\n$tbn{$word[0]}\n";        ###### The new seqID is for satisfying the requirement of ID length in genewise.  
                close DNA;
                
                open PEP,">$outdir/$count.pep" or die "$!";
                print PEP ">$word[1]\n$db{$word[1]}\n";
                close PEP;
                
                $pairs{"$outdir/$count.fa"} = "$outdir/$count.pep";
                $uniq{$word[0]} =  1;
                $count++;
            }
        }
    }
    close BLAX;
    return %pairs;
}

sub Merge_regions_result
{
    my $in = shift;
    my ($c_index,$s_index,$e_index) = ($_[0] - 1 , $_[1] - 1 , $_[2] - 1);   ### chr start end
    my (%P_chr_start,%N_chr_end);   ### split the strand -/+
    print "Do Merge process:",int(@{$in}),"\n";
    Run_jobs::Print_time();
    foreach my $line(@{$in})
    {
        my @word = split /\t/,$line;
        my @region = ($word[$c_index],$word[$s_index],$word[$e_index]);
        my $strand;
        ####Check strand
        if($region[1] <= $region[2])
        {
            $strand = "+";
        }
        else
        {
            $strand = "-";
        }
        push @region , $strand;
        
        if($region[3] eq '+')
        {
            unless($P_chr_start{$region[0]}{$region[1]})
            {
                $P_chr_start{$region[0]}{$region[1]} = (join "\t" , @region);
            }
            else
            {
                my @compar = split /\t/ , $P_chr_start{$region[0]}{$region[1]};
                &Compare_region(\@region,\@compar,'merge');
                $P_chr_start{$region[0]}{$region[1]} = (join "\t" , @compar);
            }
        }
        else
        {
            unless($N_chr_end{$region[0]}{$region[2]})
            {
                $N_chr_end{$region[0]}{$region[2]} = (join "\t", @region);
            }
            else
            {
                my @compar = split /\t/, $N_chr_end{$region[0]}{$region[2]};
                &Compare_region(\@region,\@compar,'merge');
                $N_chr_end{$region[0]}{$region[2]} = (join "\t", @region);
            }
        }
    }
#print "merging success\n";
#   Run_jobs::Print_time();
    my %output;
    my @chr = keys %P_chr_start;
    for(my $m = 0;$m < @chr;$m++)
    {
#print "Comparing $chr[$m]\n";
        my @P_out = &Compare_cycle(\%P_chr_start,$chr[$m]);
#        print "Finished $chr[$m] P_strand\n";
        my @N_out = &Compare_cycle(\%N_chr_end,$chr[$m]);
#        print "Finished $chr[$m] N_strand\n";
        @{$output{$chr[$m]}} = (@P_out,@N_out);
   }
    return %output;   #### $hash{$chr}->@regions; @regions = ($chr,$start,$end,$strand)
}

sub Compare_cycle
{
    my ($hash,$chr) = @_;
    my @sort_pos = sort {$a <=> $b}(keys %{$$hash{$chr}});
#    print "Comparing sort completed\tnumber:",int(@sort_pos),"\n";
    my @save;
    for(my $i = 0;$i < @sort_pos;$i++)
    {
        my @region = split /\t/,$$hash{$chr}{$sort_pos[$i]};
        for(my $j = $i+1;$j<@sort_pos;$j++)
        {
            my @compar = split /\t/,$$hash{$chr}{$sort_pos[$j]};
            my $judge = &Compare_region(\@compar,\@region,'merge');
            if($judge == 1)
            {
                @sort_pos = &Delete_element_from_array(\@sort_pos,$j);
#                print "$i\t$j\t",int(@sort_pos),"\n";
                $j--;
            }
        }
        push @save,(join "\t",@region);
    }
    return @save;
}

sub Delete_element_from_array
{
    my ($array,$index) = @_;
    my ($end,$start) = ($index - 1,$index + 1);
    my $num = int(@{$array}) - 1;
    my @new;
    if($index == 0)
    {
        @new = @{$array}[1..$num];
    }
    elsif($index == $num)
    {
        @new = @{$array}[0..$index-1];
    }
    else
    {
        @new = (@{$array}[0..$end],@{$array}[$start..$num]);
    }
    return @new;
}

sub Compare_region
{
    my ($new,$compar,$style) = @_; 
    $style = "compare" unless($style eq 'merge');
    my $judge = 0;
    if($$new[0] eq $$compar[0]) #Check chr;
    {
        if($$new[3] eq $$compar[3]) #Check strand;
        {
            if($$compar[3] eq '+') 
            {
                if($$compar[1] >= $$new[1] and $$compar[1] < $$new[2])
                {
                     $$compar[1] = $$new[1] if($style eq "merge");
                     $judge = 1;
                }
                if($$compar[2] > $$new[1] and $$compar[2] <= $$new[2])
                {
                    $$compar[2] = $$new[2] if($style eq "merge");
                    $judge = 1;
                }
                if($$compar[1] <= $$new[1] and $$compar[2] >= $$new[2])
                {
                    $judge = 1;
                }
            }
            else
            {
                if($$compar[1] >= $$new[2] and $$compar[1] < $$new[1])
                {
                    $$compar[1] = $$new[1] if($style eq "merge");
                    $judge = 1;
                }
                if($$compar[2] > $$new[2] and $$compar[2] <= $$new[1])
                {
                    $$compar[2] = $$new[2] if($style eq "merge");
                    $judge = 1;
                }
                if($$compar[2] <= $$new[2] and $$compar[1] >= $$new[1])
                {
                    $judge = 1;
                }
            }
        }#else{print "Check_strand\n";}
    }#else{print "Check_Chr\n";}
    return $judge;
}

sub Genewise_flow
{
    my ($pairs,$opt) = @_;
    my ($outdir,$wise_path) = ($$opt{'Outdir'},"$$opt{'Outdir'}/genewise");
    my $prefix = ($$opt{'Prefix'}) ? $$opt{'Prefix'} : "Plant";
    print "Output $prefix.*\n";
    my $ref = $$opt{'Ref'};
    open SH,">$outdir/call_genewise.sh" or die "$!";
    my @fa = keys %{$pairs};
    my @output;
    for(my $i = 0; $i < @fa; $i++)
    {
        my $num = $i+1;
        print SH "$Bin{'genewise'}   $$pairs{$fa[$i]}   $fa[$i]   -genesf  -splice_gtag -silent >  $wise_path/$num.gw\n";
        push @output,"$wise_path/$num.gw";
    }
    close SH;
#system("sh $outdir/call_genewise.sh");
    Run_jobs::Run_jobs("$outdir/call_genewise.sh",$runtype,$p_num,$queue);
    if(-e "$outdir/$prefix.gff")
    {
        system("rm $outdir/$prefix.*");
    }
##### Summary result
    my ($gff,$mut) = &Generate_genewise_summary(\@output,$prefix,$$opt{'Outdir'});
#####Output fasta
    my ($min_len,$max_ratio) = Configure::Load_opt("Fa_len","Max_N");
    my %pep = Fasta_deal::Gff_to_Pep($gff,$ref);
    my $out_fa = "$$opt{'Outdir'}/$prefix.pep";  ### the path is different from $outdir;
    Fasta_deal::Print_fasta(\%pep,$out_fa,$min_len,$max_ratio);
}

sub Generate_genewise_summary
{
    my ($file,$prefix,$outdir) = @_;
    my ($out_gff,$out_mut);
    for(my $i = 0; $i < @{$file}; $i++)
    {
        if(-e $$file[$i])
        {
            ($out_gff,$out_mut) = &Brief_genewise_result($$file[$i],$prefix,$outdir)
        }
        else
        {
            print STDERR "Warning: Can't find the result $$file[$i]\n";
        }
    }
    return ($out_gff,$out_mut);
}

sub Brief_genewise_result
{
    my ($file,$prefix,$outdir) = @_;
    my @info;   ### save location info
    my $seq_bg;
    open IN,"<$file" or die "$!";
    $/ = "//";
    chomp (my $block = <IN>);
    chomp (my $genesf = <IN>);
    close IN;
    $/ = "\n";
    ####get ID
    $block =~ /Query\s+protein:\s+(.+?)\n.+Target\s+Sequence\s+(.+?)\nStrand:\s+.+?\n/s;
    my ($pepID , $faID) = ($1 , $2);
#print "$pepID\t$faID\n";
    @info = split /_/,$faID;
    $seq_bg = $info[1];
    my $trimed_ID = (length($pepID) > 15) ? substr($pepID,0,15) : $pepID;
    ####get alignment score
    $block =~ /Score\s+(\S+)\s+bits\s+over\s+entire\s+alignment/;
    my $score = $1;
    ####Load alignment block 
    my @line = split /\n/,$block;
    my @align_block;        ##($query_pep, $match_pep, $target_pep, $nuc_line1, $nuc_line2, $nuc_line3)
    for (my $i = 0; $i < @line; $i++)
    {
        if ($line[$i] =~ /^$trimed_ID\s+/) 
        {
            for(my $j = 0; $j < 6; $j++)
            {
                substr($line[$i + $j], 0, 21) = ""; 
                $align_block[$j] .= $line[$i + $j];
            }
            $i += 6;
        }
    }
     ####Load gene structure with supporting evidence
    my @exon_location;
    my $align_len;
    my $identity;
    my @genes;
    my $gene_num = 0;
    @line = split /\n/,$genesf;
    for(my $i = 0; $i < @line; $i++)
    {
        if($line[$i] =~ /Gene\s+(\d+)$/)
        {
            $gene_num = $1;   ### index start from 1;
=cut
            if($gene_num > 1)
            {
                my $last_gene = $gene_num - 1;
                $match_pep =~ s/\s+|\+//g;
                $identity = sprintf "%.2f", length($match_pep) / $align_len * 100;
                for(my $n = 0; $n < @{$genes[$last_gene]}; $n++)
                {
                    $genes[$last_gene][$n] .= "\t$identity";
                }
            }
=cut
        }
        if($line[$i] =~ /Gene\s+(\d+)\s+(\d+)/)
        {
            $info[1] = &Trans_coord($seq_bg,$info[3],$1);
            ($info[1],$info[2]) = &Trans_region_pos($info[1],$info[3],abs($2 - $1 + 1));
            push @{$genes[$gene_num]} , "$info[0]\tGene$gene_num\tmRNA\t$info[1]\t$info[2]\t.\t$info[3]\tphase\t$faID\t";   #### Gene line #####
        }
        elsif($line[$i] =~ /\s+Exon\s+(\d+)\s+(\d+)\s+phase\s+(\d+)/)
        {
            my ($exon_bg,$exon_ed,$phase) = ($1,$2,$3);
            ($exon_bg,$exon_ed) = (&Trans_coord($seq_bg,$info[3],$exon_bg) , &Trans_coord($seq_bg,$info[3],$exon_ed));
            push @exon_location, [$exon_bg, $exon_ed, $phase];
            push @{$genes[$gene_num]} , "$info[0]\tGene$gene_num\tCDS\t$exon_bg\t$exon_ed\t.\t$info[3]\t$phase\t$faID\t";  #### part 1 each exon line####
        }
        elsif($line[$i] =~ /\s+Supporting\s+\d+\s+\d+\s+(\d+)\s+(\d+)/)
        {
            $align_len += (abs($1 - $2) + 1);
        }
    }
    #####Calculate indentity#######
    my $match_pep = $align_block[1];
    $match_pep =~ s/\s+|\+//g;
    my $match = sprintf "%.2f", length($match_pep) / $align_len * 100;
    $identity = sprintf "%.2f", length($match_pep) / length($align_block[0]);
    @align_block = &Delete_element_from_array(\@align_block,1);
    ##### Output gff.
    open OGFF,">>$outdir/$prefix.gff" or die "$!";
    for(my $i = 1 ;$i < $gene_num + 1;$i++)
    {
        for(my $j = 0; $j < @{$genes[$i]}; $j++)
        {
            print OGFF "$genes[$i][$j]\t$match\t$identity\n"; 
        }
    }
    close OGFF;
    ####Mutation type
    my @mut_sum = &Genewise_mutation_type(\@align_block,\@info,\@exon_location);
    open OMUT,">>$outdir/$prefix.mutation" or die "$!";
        for(my $j = 0;$j < @mut_sum;$j++)
        {
            print OMUT "$faID\t$mut_sum[$j]\n";
        }
    close OMUT;
    return ("$outdir/$prefix.gff","$outdir/$prefix.mutation");
}

sub Genewise_mutation_type
{
    my ($align,$nuc_info,$exon) = @_;
    my $length = length($$align[0]);
    my @mut_output;
    my ($bg,$strand) = ($$nuc_info[1],$$nuc_info[3]);
    my $exon_count = 0;
    for(my $i = 0; $i < $length; $i++)
    {
        my @char;       #(qurey_pep,target_pep,$nc_line1,$nc_line2,$nc_line3)
        my ($nc_bg,$nc_ed,$mut_char,$mut_type);
        for(my $j = 0; $j < @{$align}; $j++)
        {
            $char[$j] = substr($$align[$j],$i,1)
        }
        if($char[0] eq '-')
        {
            ($nc_bg,$nc_ed) = ($char[2] =~ /\d+/) ? &Trans_region_pos($bg,$strand,$char[2]) : &Trans_region_pos($bg,$strand,3);
            $mut_char = join "",(@char[2..4]);
            $mut_type = 'I';
        }
        elsif($char[1] eq '-')
        {
            ($nc_bg,$nc_ed) = ($strand eq '+') ? &Trans_region_pos(($bg-1),$strand,2):&Trans_region_pos(($bg+1),$strand,2);
            $mut_char = $char[0];
            $mut_type = 'D';
        }
        elsif($char[1] eq 'X' && $char[0] ne 'U' && $char[0] ne 'X')
        {
            ($nc_bg,$nc_ed) = &Trans_region_pos($bg,$strand,3);
            my $codon = join "",(@char[2..4]);
            $mut_char = &Check_codon_mutation($codon,$char[0]);
            $mut_type = 'S';
        }
        elsif($char[1] eq '!')
        {
            die "first line of nuc isn't numberic : $char[2]\n" if($char[2] !~ /\d+/);
            ($nc_bg,$nc_ed) = &Trans_region_pos($bg,$strand,$char[2]);
            $mut_char = "$char[2]:$char[0]";
            $mut_type = "F";
        }
        else
        {
            #print (join "\n",@char);
        }
        if($mut_type)
        {
            push @mut_output, "$nc_bg\t$nc_ed\t$mut_char\t$mut_type";
        }
        #####Move the position of nuc
        if($char[1] =~ /\s+/ && $char[3] eq "<")
        {
            $i += 22;
            $exon_count++;
            $bg = $$exon[$exon_count] -> [0];
        }
        elsif($char[1] eq '!')
        {
#           my $next_ba = substr($$align[2], $i+1, 1);
#           my $next_ba2 = substr($$align[3], $i+1, 1);
#           $exon_count++ unless( ($next_ba=~/[ATCGN]/ && $next_ba2!~/[atcgnN]/) || $next_ba =~ /\d+/ );
            if ($char[2] =~ /\d+/) 
            {
                $bg = &Trans_coord($bg,$strand,$char[2]);
            }
            else
            {
#                $bg = $$exon[$exon_count]->[0];
                next;
            }
        }
        elsif($char[1] eq '-')
        {
            next;
        }
        else
        {
            if($char[2] =~ /[atcgn]/)
            {
                $bg = &Trans_coord($bg,$strand,3);
            }
            elsif($char[2] =~ /N/ && $char[3] =~ /[atcgnN]/ && $char[4] =~ /[atcgnN]/)
            {
                $bg = &Trans_coord($bg,$strand,3);
            }
            elsif($char[2] =~ /[ATGCN]/)
            {
                $bg = &Trans_coord($bg,$strand,1);
            }
        }
    }
    return @mut_output;
}

sub Check_codon_mutation
{
    my ($codon,$aa) = @_;
    $codon = uc($codon);
    my @ref = Fasta_deal::AA_to_Nucl($aa);
    my %mut;
    my %mut_num;
    for(my $n = 0;$n < @ref;$n++)
    {
        for(my $i = 0;$i < 3;$i++)
        {
            my $q_base = substr($codon,$i,1);
            my $r_base = substr($ref[$n],$i,1);
            my $site = $i + 1;
            if($q_base ne $r_base)
            {
                push @{$mut{"$ref[$n]>$codon"}},"$site:$r_base>$q_base";
                $mut_num{"$ref[$n]>$codon"}++;
            }
        }
    }
    foreach my $r_codon (sort {$mut_num{$a} <=> $mut_num{$b}} keys %mut_num)
    {
        my $info = join ",",@{$mut{$r_codon}};
        return "$r_codon:$info";
        last;
    }
}
sub Trans_region_pos
{
    my ($start,$strand,$length) = @_;
    my $end = &Trans_coord($start,$strand,$length);
    return ($start,$end);
}

sub Trans_coord
{
    my($pos,$strand,$length) = @_;
    my $out;
    if($strand eq '+')
    {
        $out = $pos + $length - 1;
    }
    elsif($strand eq '-')
    {
        $out = $pos - $length + 1;
    }
    else
    {
        die "Can't recognise the strand info : $strand\n";
    }
    return $out;
}

sub Check_with_gff
{
    my ($region,$id_file,$gff) = @_;   ###$region->%hash;$id->file_address;$gff->file_address.
    ####read input for gene_id
    my @id;
    open ID,"<$id_file" or die "$!";
    while(<ID>)
    {
        next unless($_=~s/^>//);
        chomp;
        my @word = split /\s+/,$_;
        push @id,$word[0];
    }
    close ID;
    ####search gff to locate gene.
    if($gff =~/\.gz$/)
    {
        open GFF,"gzip -dc $gff |" or die "$!";
    }
    else
    {
        open GFF,"<$gff" or die "$!";
    }
    my $max_overlap = Configure::Load_opt('Overlap');
    while(<GFF>)
    {
        next if($_=~/^#/);
        my @word = split /\t/;
        my $line = $_;
        next unless($word[2]=~/mRNA/i or $word[2]=~/gene/i);
        for(my $i = 0;$i < @id;$i++)
        {
            if($line =~/$id[$i]/i)
            {
                my @gff_region;
                if($word[6] eq '+')
                {
                    @gff_region = ($word[0],$word[3],$word[4],$word[6])
                }
                elsif($word[6] eq '-')
                {
                    @gff_region = ($word[0],$word[4],$word[3],$word[6])
                }
                else
                {
                    die "Check_with_gff(module_name)::Dose The 7th column mean the strands?\n";
                }
                if(!$$region{$word[0]})
                {
                    print "$word[0]\n";
                    next;
                }
                for(my $j = 0;$j < @{$$region{$word[0]}};$j++)
                {
                    my @pos = split /\t/,${$$region{$word[0]}}[$j];
                    my $judge = &Compare_region(\@pos,\@gff_region,'merge');   
                    if($judge == 1)
                    {
                        my $gene_len = abs($word[4] -$word[3]) + 1;
                        my $new_len = abs($gff_region[2] - $gff_region[1]) + 1;
                        my $pos_len = abs($pos[2] - $pos[1]) + 1;
                        my $overlap_len = ($pos_len - ($new_len - $gene_len));
                        my $overlap_ratio = $overlap_len / $pos_len;
                        if($overlap_ratio > $max_overlap || $overlap_len == $pos_len)
                        {
#                            print "Check:$word[3]\t$word[4]\t$pos[1]\t$pos[2]\t$overlap_len\t$pos_len\t$overlap_ratio\t$max_overlap\n";
                            @{$$region{$word[0]}} = &Delete_element_from_array(\@{$$region{$word[0]}},$j);
                            $j--;
                        }
                    }
                }
            }
        }      
    }
    close GFF;
    my @chr =keys %{$region};
    my $count = 0;
    foreach(@chr)
    {
        $count += int(@{$$region{$_}});
    }
    print "Second_gff_check confirm:$count\n";
    Run_jobs::Print_time();
}
                
    


    
    

1;#keep require happy
