#!/usr/bin/perl -w

package Configure;
use strict;
use FindBin qw($Bin);
use Exporter;
our @EXPORT = qw(Bin_Adr Load_opt);
our @EXPORT_OK  = qw(%Bin);

our %Bin = (
        'formatdb' => "",
        'hmmsearch' => "",
        'blastp' => "",
        'blastx' => "",
        'tblastn' => "",
        'genewise' => "",
        'hmmbuild' => "",
        'clustalw2' => "",
        );
our %opt_save;

sub Save_opt
{
    my $opt = shift;
    %opt_save = %{$opt};
}

sub Load_opt
{
    my @name = @_;
    my @output;
    for(my $i = 0; $i < @name; $i++)
    {
        if($opt_save{$name[$i]})
        {
            $output[$i] = $opt_save{$name[$i]};
        }
        else
        {
            if($name[$i] eq 'Process')
            {    
                $output[$i] = 1;
            }
            elsif($name[$i] eq 'Retain')
            {
                $output[$i] = 'F';
            }
            elsif($name[$i] eq 'AA_len')
            {
                $output[$i] = '';
            }
            elsif($name[$i] eq 'Identity')
            {
                $output[$i] = '';
            }
            elsif($name[$i] eq 'Max_gap')
            {
                $output[$i] = '';
            }
        }
    }
    if(int(@output) == 1)
    {
        return $output[0];
    }
    else
    {
        return @output;
    }
}

sub Check_Configure
{
    if(-e "$Bin/Configure.txt")
    {
        &Load_Configure(\%Bin,"$Bin/Configure.txt");
        my $judge = &Check_missing(\%Bin);
        if($judge)
        {
            &Creat_file(\%Bin,"$Bin/Configure.txt");
            die "\nSome tools are required for this pipeline\nPlease fill the information in $Bin/Configure.txt\n\n\n";
        }
    }
    else
    {
        &Creat_file(\%Bin,"$Bin/Configure.txt");
        die "\nSome tools are required for this pipeline\nPlease fill the information in $Bin/Configure.txt\n\n\n";
    }
}

sub Bin_Adr
{
    &Check_Configure();
    return %Bin;
}

sub Load_Configure
{
    my ($bin,$in) = @_;
    open CONF,"<$in" or die "$!";
    while(<CONF>)
    {
        next if($_=~/^#/);
        chomp;
        $_ =~s/\s+//g;
        my @info = split /=/,$_;
        if(exists $$bin{$info[0]})
        {
            $$bin{$info[0]} = $info[1];
        }
        else
        {
            my $wanted_list = join "\n",(keys %{$bin});
            die "Maybe the $info[0] is a typo\n Please check in $Bin/Configure.txt\nThe wanted list:\n$wanted_list\n";
        }
    }
    close CONF;
}

sub Creat_file
{
    my ($bin,$out) = @_;
    my @wanted = keys %{$bin};
    open OUT,">$out" or die "$!";
    for(my $i = 0;$i < @wanted;$i++)
    {
        if($$bin{$wanted[$i]})
        {
            print OUT "$wanted[$i]=$$bin{$wanted[$i]}\n";
        }
        else
        {
            print OUT "$wanted[$i]=\n";
        }
    }
    close OUT;
}

sub Check_missing
{
    my $bin = shift;
    my @program = keys %{$bin};
    my $recreate = 0;
    for(my $i = 0;$i < @program;$i++)
    {
        unless($$bin{$program[$i]})
        {
            print STDERR "The program : #$program[$i]# is not available\n";
            $recreate = 1;
        }
    }
    return $recreate;
}

1;
