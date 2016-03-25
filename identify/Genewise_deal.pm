#!usr/bin/perl -w
use strict;

package Genewise_deal;

sub Genewise_flow
{
    my ($pairs,$loca,$opt) = @_;
    my ($outdir) = "$$opt{'Outdir'}/genewise";
    open SH,">$outdir/call_genewise.sh" or die "$!";
    my @fa = keys %{$pairs};
    my @output;
    for(my $i = 0; $i < @fa; $i++)
    {
        my $num = $i+1;
        print SH "$Bin{'genewise'}   $$pairs{$fa[$i]}   $fa[$i]   -genesf  -splice_gtag  >  $outdir/$num.gw\n";
        push @output,"$outdir/$num.gw";
    }
    close SH;
    Run_jobs::Run_jobs("$outdir/call_genewise.sh",$opt);
    &Generate_genewise_summary(\@output,$loca,$opt);
}

sub Generate_genewise_summary
{
    my ($file,$loca,$opt) = @_;
    my $outdir = "$$opt{'Outdir'}/genewise";
    my $prefix = ($$opt{'Prefix'}) ? $$opt{'Prefix'} : "Plant";
    for(my $i = 0; $i < @{$file}; $i++)
    {
        if(-e $$file[$i])
        {
            my ($alg,$mut,$gff) = &Brief_genewise_result($$file[$i],$$loca[$i])
        }
        else
        {
            print STDERR "Warning: Can't find the result $$file[$i]\n";
        }
    }
}

sub Brief_genewise_result
{
    my ($file,$loca) = @_;
    my @info = split /_/,$loca;
    my $seq_bg = $info[1];
    open IN,"<$file" or die "$!";
    $/ = "//";
    while(<IN>)
    {
        chomp;
        ####get ID
        $_ =~ /Query\s+protein:\s+(.+?)\n.+Target\s+Sequence\s+(.+?)\nStrand:\s+.+?\n/s;
        my ($pepID , $fqID) = ($1 , $2);
        my $trimed_ID = (length($pepID) > 15) ? substr($pepID,0,15) : $pepID;
        ####get output aligned block
        $_ =~ /.+genewise\s+output(.+)\/\//s;
        my $block = $1;
        ####get alignment score
        $block =~ /Score\s+(\S+)\s+bits\s+over\s+entire\s+alignment/;
        my $score = $1;
        ####Load alignment block 
        my @line = split /\n/,$block;
        my @align_block;        ##($query_pep, $match_pep, $target_pep, $nuc_line1, $nuc_line2, $nuc_line3)
        for (my $i = 0; $i < @line; $i++)
        {
            if ($line[$i] =~ /^$trimed_ID/) 
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
        my $genesf = <IN>;
        chomp $genesf;
        my @nuc_align_block;
        my @pep_align_block;
        my @exon_location;
        my $match_len;
        my $identity;
        @line = split /\n/,$genesf;
        for(my $i = 0; $i < @line; $i++)
        {
            if($line[$i] =~ /Gene\s+(\d+)\s+(\d+)/)
            {
                $info[1] = &Trans_coord($seq_bg,$info[3],$1);
               ($info[1],$info[2]) = &Trans_region_pos($info[1],$info[3],abs($2 - $1 + 1));
            }
            elsif($line[$i] =~ /\s+Exon\s+(\d+)\s+(\d+)\s+phase\s+(\d+)/)
            {
                my ($exon_bg,$exon_ed,$phase) = ($1,$2,$3);
                ($exon_bg,$exon_ed) = (&Trans_coord($seq_bg,$info[3],$exon_bg) , &Trans_coord($seq_bg,$info[3],$exon_ed));
                push @nuc_align_block, ($exon_bg,$exon_ed);
                push @exon_location, [$exon_bg, $exon_ed, $phase];
            }
            elsif($line[$i] =~ /\s+Supporting\s+\d+\s+\d+\s+(\d+)\s+(\d+)/)
            {
                push @pep_align_block, ($1, $2);
                $match_len += (abs($1 - $2) + 1);
                my $match_pep = $align_block[1];
                $match_pep =~ s/\s+|\+//g;
                $identity = sprintf "%.2f", length($match_pep) / $match_len * 100;
            }
        }
        ####Mutation type
        @align_block = &Delete_element_from_array(\@align_block,1);  ###Delete the match_pep. useless in detecting mut type.
        my @mut_sum = &Genewise_mutation_type(\@align_block,\@info,\@exon_location);
    } ### end while()
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
            $char[$i] = substr($$align[$i],$i,1)
        }
        if($char[0] eq '-')
        {
            ($nc_bg,$nc_ed) = ($char[2] =~ /\d+/) ? &Trans_region_pos($bg,$strand,$char[2]) : &Trans_region_pos($bg,$strand,3);
            $mut_char = join "",(@char[2..4]);
            $mut_type = 'I';
        }
        elsif($char[1] eq '-')
        {
            ($nc_bg,$nc_ed) = ($$nuc_info[3] eq '+') ? &Trans_region_pos(($bg - 1),$strand,2) : &Trans_region_pos(($bg + 1),$strand,2);
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
        push @mut_output, "$nc_bg\t$nc_ed\t$mut_char\t$mut_type";
        #####Move the position of nuc
        if($char[1] =~ /\s+/ && $char[3] eq "<")
        {
            $i += 22;
            $exon_count++;
            $bg = $$exon[$exon_count] -> [0];
        }
        elsif($char[1] eq '!')
        {
            my $next_ba = substr($$align[2], $i+1, 1);
            my $next_ba2 = substr($$align[3], $i+1, 1);
            $exon_count++ unless( ($next_ba=~/[ATCGN]/ && $next_ba2!~/[atcgnN]/) || $next_ba =~ /\d+/ );
            if ($next_ba =~ /\d+/) 
            {
                $bg = &Trans_coord($bg,$strand,$char[2]);
            }
            else
            {
                $bg = $$exon[$exon_count]->[0];
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

