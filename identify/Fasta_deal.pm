package Fasta_deal;
use strict;
use warnings;
use File::Basename;
use File::Path; ##mkpath() 

sub Print_fasta
{
    my ($fa,$out,$min_len,$max_ratio) = @_;
    $min_len ||= 0;  ####default config will disable filtration. 
    $max_ratio ||= 1;
    open FAOUT,">$out";
    my @id = keys %{$fa};
    my $before_filter = int(@id);
    my $after_filter;
    for(my $i = 0;$i < @id;$i++)
    {
        my $seq = $$fa{$id[$i]};
        my $len = length ($seq);
        my $n += $seq=~s/[nN\-\?\*]//g;
        my $n_ratio = $n / $len;
        if($len > $min_len && $n_ratio < $max_ratio)
        {
            print FAOUT ">$id[$i]\n$$fa{$id[$i]}\n";
            $after_filter++;
        }
    }
    print STDERR "Filter_criterion: min_len:$min_len\tN_ratio:$max_ratio\nBefore:$before_filter\tAfter:$after_filter\n";
    close FAOUT;
}

        

sub Gff_to_Pep
{
    my ($gff,$ref) = @_;

    if($gff =~/\.gz$/)
    {
        open GFF,"gzip -dc $gff|" or die "$!";
    }
    else
    {
        open GFF,"<$gff" or die "$!";
    }
    my %seq = &Load_fasta_by_id($ref);
    my %cds;
    my $ID;
    while(<GFF>)
    {
        chomp;
        if($_ =~ /^#/)
        {
            next;
        }
        my @col = split /\t/,$_;
        if($col[2] =~ /mRNA/i)
        {
            my @region = ($col[0],$col[3],$col[4],$col[6]);
            $col[-1] =~ /ID=(.+?);/;
            $ID = $1;
            if(!$ID)
            {
                $ID = join "_",@region;
                $ID .="_$col[1]"; 
#die "Can't get ID from mRNA line. /ID=(.+?)/. Check formate:\t$col[-1]\n";
            }
        }
        elsif($col[2] =~ /CDS/i)
        {
            my @region = ($col[0],$col[3],$col[4],$col[6]);
            $cds{$ID} .= &Extract_subseq(\%seq,@region);
        }
    }
    my %pep = &Trans_nuc_to_aa(\%cds);
    return %pep;
}

sub Trans_nuc_to_aa
{
    my $cds = shift;
    my @id = keys %{$cds};
    my %trans = &Reverse_ref_codon();
    my %pep;
    my %error;
    for(my $i = 0;$i < @id;$i++)
    {
        my $seq = $$cds{$id[$i]};
        my $len = length($seq);
        next if($len <= 102);
#&Check_strcuture($seq);
        if(($len % 3) > 0)
        {
            die "the cds length isn't triple\n";
        }
        else
        {
            for(my $n = 0; $n < $len; $n += 3)
            {
                my $codon = substr($seq,$n,3);
                ####Check Structure
                if($n == 0)
                {
                    if($trans{$codon} ne "M")
                    {
                        print STDERR "$id[$i]:Lack start\n";
                    }
                }
                if($n == ($len - 3))
                {
                    if($trans{$codon} ne "*")
                    {
                        print STDERR "$id[$i]:Lack stop\n";
                    }
                }
                if($trans{$codon} and $trans{$codon} ne "*")
                {
                    $pep{$id[$i]} .= $trans{$codon}
                }
                elsif($trans{$codon} and $trans{$codon} eq "*")
                {
                    if($n != ($len - 3))
                    {
                        print STDERR "$id[$i]:Inside$n:Stop\n";
                    }
                }
                else
                {
                    print STDERR "$codon can't be recognise\n";
                }
            }
        }
    }
    return %pep;
}

sub Check_structure
{
    my %trans = &Reverse_ref_codon();
    my $seq = shift;
    my $start = substr($seq,0,3);
    my $last = length($seq) - 1; ##count start at 0;
    my $end = substr($seq,($last - 3),3);
    my @output;
    if($trans{$start} eq 'M')
    {
        $output[0] = 1;
    }
    else
    {
        $output[0] = 0;
    }
    if($trans{$end} eq '*')
    {
        $output[1] = 1;
    }
    else
    {
        $output[1] = 0;
    }
    return @output;
}


sub Load_gff_byID
{
    my ($gff,$fa) = @_;
    my %fasta = &Load_fasta_by_id($fa);
    my @id = keys %fasta; undef %fasta;
    my %output;
    if($gff =~/\.gz$/)
    {
         open GFF,"gzip -dc $gff |" or die "$!";
    }
    else
    {
         open GFF,"<$gff" or die "$!";
    }
    while(<GFF>)
    {
        next if($_=~/^#/);
        my @word = split /\t/;
        next unless($word[2]=~/mRNA/i or $word[2]=~/gene/i);
        my $line = $_;
        for(my $i = 0;$i < @id;$i++)
        {
             if($line =~/$id[$i]/i)
             {
                 my $gff_region;
                 if($word[6] eq '+')
                 {
                      $gff_region = join "\t",($word[0],$word[3],$word[4],$word[6]);
                 }
                 elsif($word[6] eq '-')
                 {
                     $gff_region = join "\t",($word[0],$word[4],$word[3],$word[6]);
                 }
                 else
                 {
                     die "Dose The 7th column of Gff is the strands?\n";
                 }
                 push @{$output{$word[0]}},$gff_region;
             }
        }
    }
    return %output;
}


sub Extract_by_regions    
{
    my ($region,$ref,$out) = @_;    ###the $region is a pointer specifically point to the data structure : @{$hash{chr}};
    my %ref = &Load_fasta_by_id($ref);
    open OUT,">$out" or die "$!";
    my @chr_id = keys %{$region};
    for(my $i = 0;$i < @chr_id;$i++)
    {
        for(my $j = 0; $j < @{$$region{$chr_id[$i]}};$j++)
        {
            my @region = split /\t/,$$region{$chr_id[$i]}->[$j];
            my $sub_seq = &Extract_subseq(\%ref,@region);
            my $name = join "_" , @region;
            print OUT ">$name\n";
            print OUT "$sub_seq\n";
        }
    }
    close OUT;
}

sub Load_fasta_by_id
{
    my $in = shift;
    $/ = ">";
    my %save;
    if($in =~/\.gz$/)
    {
        open IN,"gzip -dc $in | " or die "$!";
    }
    else
    {
        open IN,"<$in" or die "$!";
    }
    <IN>; #First lines will be ">";
    while(<IN>)
    {
        chomp;
        my @lines = split /\n/;
        my $id = shift @lines;
        if($id =~/\s/)
        {
            my @word = split /\s+/,$id;
            $id = $word[0];
        }
        $save{$id} = (join "",@lines);
    }
    close IN;
    $/ = "\n";
    return %save;
}

sub Extract_subseq
{
    my ($ref,$chr,$start,$end,$strand) = @_;
    $strand ||= "+";
    die "the chr id don't match to the ref : $chr\nFrom &Extract_subseq\n" if(!$$ref{$chr});
    my $sub_seq;
    if($strand eq "+"){
        $sub_seq = substr($$ref{$chr},$start-1,$end-$start+1);
    }
    elsif($strand eq "-")
    {
        $sub_seq = substr($$ref{$chr},$end-1,$start-$end+1);
        $sub_seq = &Trans_strand($sub_seq);
    }
    else
    {
        die "the strand must be +/-\n";
    }
}

sub Trans_strand
{
    my %trans=(
            "A" => "T",
            "T" => "A",
            "C" => "G",
            "G" => "C",
            "N" => "N",
            "a" => "t",
            "t" => "a",
            "c" => "g",
            "g" => "c",
            );
    my $in = shift;
    my @base = split //,$in;
    @base = &Random_pick_het(\@base);
    @base = map{$trans{$_}} @base;
    @base = reverse(@base);
    my $trans_seq = join "",@base;
    return $trans_seq;
}

sub Random_pick_het
{
	my $geno = shift;
	my %het = &Het_code();
	my @output;
	my $judge = 0;
	for(my $i = 0; $i < @{$geno}; $i++)
	{
		my $trans;
		if($het{$$geno[$i]})
		{
			my $num = int(@{$het{$$geno[$i]}});
			my $random = rand($num);
			my $index = int($random);
			$trans = ${$het{$$geno[$i]}}[$index];
			$judge = 1;
		}
		else
		{
			$trans = $$geno[$i];
		}
		$output[$i] = $trans;
	}
	print "Warning:Seq have Het code\n" if($judge);
	return @output;
}
					
sub Split_fasta_average_size
{
    my ($Cutf,$in_fa,$Outdir) = @_;
	my $Total_len=0;
    if($in_fa =~/\.gz$/)
    {
        open IN,"gzip -dc $in_fa |" or die "$!";
    }
    else
    {
        open IN,"<$in_fa" or die "$!";
	}
    $/=">";<IN>;$/="\n";
	while(<IN>) 
    {
		chomp;
		my $head=$_;
		$/=">";
		my $seq=<IN>;
		chomp $seq;
		$/="\n";
		$seq=~s/\s//g;
		$Total_len+=length($seq);
    }
	close IN;
	##warn $Total_len;
	my $file_name = basename $in_fa; ## only file name of $infile
	chomp $file_name;
	my $out_dir = (defined $Outdir) ? "$Outdir/$file_name.cut" : "./$file_name.cut";
	if(-e "$out_dir")
    {
        system("rm  $out_dir/*.fa");
    }
    else
    {
        mkpath($out_dir);
    }
	my $Sub_len = int($Total_len / $Cutf);
	##warn $Sub_len;
	my $Cur_len = 0;
	my $Cur_content;
	my $file_mark = 1;
    if($in_fa =~/\.gz$/)
    {
        open IN,"gzip -dc $in_fa |" or die "$!";
    }
    else
    {
        open IN,"<$in_fa" or die "$!";
    }
	$/ = ">";<IN>;$/ = "\n";
    my @subfiles;   ###return values
	while(<IN>) 
    {
		chomp;
		my $head=$_;
		$/=">";
		my $seq=<IN>;
		chomp $seq;
		$/="\n";
		my $str = $seq;
		$str=~s/\s//g;
		$Cur_len+=length($str);
		$Cur_content .= ">$head\n$seq";
		if ($Cur_len >= $Sub_len) 
        {
			open OUT,">$out_dir/$file_name.$file_mark.fa" || die "fail";
			push @subfiles,"$out_dir/$file_name.$file_mark.fa";
            print OUT $Cur_content;
			close OUT;
			$Cur_content = "";
			$Cur_len = 0;
			$file_mark++;
		}
    }
    close IN;
##make the last file;
	if($Cur_content) 
    {
	    open OUT,">$out_dir/$file_name.$file_mark.fa" || die "fail";
	    push @subfiles,"$out_dir/$file_name.$file_mark.fa";
        print OUT $Cur_content;
	    close OUT;
	}
    return @subfiles;
}

sub Extract_fa_byID
{
	my ($id,$in,$out) = @_;
	open IN,"<$in" or die "$!";
	open OUT,">$out" or die "$!";
	$/=">";
	<IN>;
	while(<IN>)
	{
		chomp;
		my @line = split /\n/;
		my $name = shift @line;
		if(exists $$id{$name} && $$id{$name} == 1)
		{
			my $seq = join "\n",@line;
			print OUT ">$name\n";
			print OUT "$seq\n";
		}
	}
	$/="\n";
	close IN;
	close OUT;
}

sub AA_to_Nucl
{
    my $aa = shift;
    my %ref_codon = &Ref_codon();
    if($ref_codon{$aa})
    {
        return @{$ref_codon{$aa}};
    }
    else
    {
        die "Wrong Character of Acid amio : $aa\n";
    }
}

sub Reverse_ref_codon
{
    my %aa_codon = &Ref_codon();
    my @aa = keys %aa_codon;
    my %reverse;
    foreach my $temp(@aa)
    {
        foreach(@{$aa_codon{$temp}})
        {
            $reverse{$_} = $temp;
        }
    }
    return %reverse;
}


sub Ref_codon
{
    my %ref_codon;
    $ref_codon{'A'} = ['GCT','GCC','GCA','GCG'];
    $ref_codon{'R'} = ['CGT','CGC','CGA','CGG','AGA','AGG'];
    $ref_codon{'N'} = ['AAT','AAC'];
    $ref_codon{'D'} = ['GAT','GAC'];
    $ref_codon{'C'} = ['TGT','TGC'];
    $ref_codon{'Q'} = ['CAA','CAG'];
    $ref_codon{'E'} = ['GAA','GAG'];
    $ref_codon{'G'} = ['GGT','GGC','GGA','GGG'];
    $ref_codon{'H'} = ['CAT','CAC'];
    $ref_codon{'I'} = ['ATT','ATC','ATA'];
    $ref_codon{'L'} = ['CTT','CTC','CTA','CTG','TTA','TTG'];
    $ref_codon{'K'} = ['AAA','AAG'];
    $ref_codon{'M'} = ['ATG'];
    $ref_codon{'F'} = ['TTT','TTC'];
    $ref_codon{'P'} = ['CCT','CCC','CCA','CCG'];
    $ref_codon{'S'} = ['TCT','TCC','TCA','TCG','AGT','AGC'];
    $ref_codon{'T'} = ['ACT','ACC','ACA','ACG'];
    $ref_codon{'W'} = ['TGG'];
    $ref_codon{'Y'} = ['TAT','TAC'];
    $ref_codon{'V'} = ['GTT','GTC','GTA','GTG'];
    $ref_codon{'*'} = ['TAA','TAG','TGA'];
    return %ref_codon;
}

sub Het_code
{
	my %het_code;
	$het_code{'M'} = ['A','C'];
	$het_code{'R'} = ['A','G'];
	$het_code{'Y'} = ['C','T'];
	$het_code{'S'} = ['G','C'];
	$het_code{'W'} = ['A','T'];
	$het_code{'K'} = ['G','T'];
	$het_code{'B'} = ['C','G','T'];
	$het_code{'D'} = ['A','T','G'];
	$het_code{'H'} = ['A','C','T'];
	$het_code{'V'} = ['A','C','G'];
	return %het_code;
}

1;
