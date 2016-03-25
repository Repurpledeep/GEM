#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin";
push @INC,$Bin;
use identify::Identify_main;
use identify::Formate_main;

if(@ARGV < 1)
{
	&Usage_main();
}
else
{
	my $Main = $ARGV[0];
	if($Main eq "Identify")
	{
        Identify_main::Main(@ARGV);
	}
	elsif($Main eq "Formate")
    {
        Formate_main::Main(@ARGV);
    }
    else
	{
		&Usage_main();
	}
}


sub Usage_main
{
	my $output =<<Usage_main;

Program: Genfat (Gene family toolkit)
Version:0.01    repurpledeep\@gmail.com Update_date:    Step 10 2015

    Usage:

        Identify    Toolkit for identify Resistance genes.
        Formate     Deal with Formate.
Usage_main
	print "$output\n";
}
