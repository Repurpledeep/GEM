package Domains_summary;
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';

sub Usage
{
    my $usage =<<Usage;
                
                Usage:  DOMAIN -Pfam <*.Pfamscan> -Marcoil <*.CC> -Interpro <*.Interscan> -Output <*>
                        
                        -Pfam       Result from PfamScan.pl
                        
                        -Marcoil    Result from Marcoil.

                        -Interpro   Result from InterproScan.

                        -Output     Output file.

                Options:

                        -Pvalue     Threshold for PfamScan.[0.05]
                        -Cthres     Threshold for Marcoil.[50]
Usage
    print "$usage\n";
}
                    

sub Load_opt
{
    my %options;
    GetOptions(
            "Pfam:s" => \$options{'Pfam'},
            "Marcoil:s" => \$options{'Marcoil'},
            "Output:s" => \$options{'Output'},
            "Interpro:s" => \$options{'Interpro'},
            "Pvalue:i" => \$options{'Pvalue'},
            "Cthres:i" => \$options{'Cthres'},
            );
    if($options{'Pfam'} || $options{'Marcoil'} || $options{'Interpro'})
    {
        &Check_default_opt(\%options);
        &Function_summary(\%options);
    }
    else
    {
        &Usage();
    }
}

sub Check_default_opt
{
    my $opt = shift;
    $$opt{'Pfam'} = abs_path($$opt{'Pfam'}) if($$opt{'Pfam'});
    $$opt{'Marcoil'} = abs_path($$opt{'Marcoil'}) if($$opt{'Marcoil'});
    $$opt{'Interpro'} = abs_path($$opt{'Interpro'}) if($$opt{'Interpro'});
    $$opt{'Pvalue'} ||= 0.05;
    $$opt{'Cthres'} ||= 50;
}

sub Function_summary
{
    my $opt = shift;
    my %output;
    if($$opt{'Pfam'})
    {
        &Load_Pfamscan($$opt{'Pfam'},\%output);
    }
    if($$opt{'Marcoil'})
    {
        &Load_Marcoil($$opt{'Marcoil'},\%output);
    }
    if($$opt{'Interpro'})
    {
        &Load_Interproscan($$opt{'Interpro'},\%output);
    }
    if($$opt{'Output'})
    {
        &Print_to_file(\%output,$$opt{'Output'});
    }
    else
    {
        &Print_to_screen(\%output);
    }
}

sub Print_to_file
{
    my ($data , $out) = @_;
    my @id = keys %{$data};
    my %summary;
    my %output;
    for(my $i = 0; $i < @id; $i++)
    {
        my @domains = sort (keys %{$$data{$id[$i]}});
        my $name = join ":",@domains;
        push @{$output{$name}},$id[$i];
        $summary{$name}++;
    }
    open OUT,">$out" or die "$!";
    my @group = sort {$summary{$b} <=> $summary{$a}} (keys %summary);
    my @count;
    for(my $j = 0; $j < @group; $j++)
    {
        my @new_id = sort (@{$output{$group[$j]}});
        for(my $m = 0; $m < @new_id;$m++)
        {
            print OUT "$new_id[$m]\t$group[$j]\n";
        }
        push @count,$summary{$group[$j]};
    }
    my $line = join "\t",@group;
    my $line2 = join "\t",@count;
    print OUT "$line\n";
    print OUT "$line2\n";
    close OUT;
}

sub Print_to_screen
{
    my ($data , $out) = @_;
    my @id = sort(keys %{$data});
    my %summary;
    my %output;
    for(my $i = 0; $i < @id; $i++)
    {
        my @domains = sort (keys %{$$data{$id[$i]}});
        my $name = join ":",@domains;
        $summary{$name}++;
        push @{$output{$name}},$id[$i];
    }
    my @group = sort {$summary{$b} <=> $summary{$a}} (keys %summary);
    my @count;
    for(my $j = 0; $j < @group; $j++)
    {
        my @new = sort @{$output{$group[$j]}};
        for(my $n = 0;$n < @new; $n++)
        {
            print "$new[$n]\t$group[$j]\n";
        }
        push @count,$summary{$group[$j]};
    }
    my $line = join "\t",@group;
    my $line2 = join "\t",@count;
    print "$line\n";
    print "$line2\n";
}

sub Load_Pfamscan
{
    my ($in , $output , $threshold) = @_;
    $threshold ||= 1;
    open PFAM,"<$in" or die "$!";
    while(<PFAM>)
    {
        chomp;
        next if($_=~/^#/ || !($_));
        my @col = split /\s+/;
        $col[6] =~s/_\d+$//;
        if($col[12] < $threshold)
        {
            $$output{$col[0]}{$col[6]} .= "$col[3]\t$col[4]\n";
        }
    }
    close PFAM;
}

sub Load_Marcoil
{
    my ($in , $output , $threshold) = @_;
    $threshold ||= 50;
    open MARC,"<$in" or die "$!";
    my ($id,$number);
    while(<MARC>)
    {
        chomp;
        next if(!$_);
        if($_=~/^>/)
        {
            my @word = split /\s+/;
            $id = shift @word;
            $id =~s/^>//;
        }
        if($_=~/NUMBER PREDICTED COILED-COIL DOMAINS WITH THRESHOLD 2\.0 : (\d+)/)
        {
            $number = $1;
            for(my $i = 0 ;$i < $number; $i++)
            {
                chomp (my $line = <MARC>);
                $line =~/from (\d+) to (\d+) \(length = \d+\) with max = (.+)/;
                my ($start , $end , $value) = ($1, $2 , $3);
                if($value >= $threshold)
                {
                    $$output{$id}{'CC'} .= "$start\t$end\n";
                }
            }
        }
    }
    close MARC;
}
            
sub Load_Interproscan
{
    my ($in , $output , $threshold) = @_;
    open INTE,"<$in" or die "$!";
    while(<INTE>)
    {
    }
}

1; #keep require happy
