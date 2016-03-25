package Identify_main;
use strict;
use warnings;
use identify::Protein_identify;
use identify::Genome_identify;
use identify::Domains_summary;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(Main);


sub Main
{
	if(int(@_) < 2)
	{
		&Usage();
	}
	else
	{
		shift @_;
		my $Sub = shift @_;
		if($Sub eq "PROTEIN")
		{
			Protein_identify::Main(@_);
		}
		elsif($Sub eq "GENOME")
		{
			Genome_identify::Main(@_);
		}
		elsif($Sub eq "DOMAIN")
        {
            Domains_summary::Load_opt(@_);
        }
        else
		{
			&Usage();
		}
	}
}



sub Usage
{
    my $output =<<USAGE;

        Identify Usage:
                PROTEIN     Identify protein from pep file.
                GENOME      Identify protein from genome.
                DOMAIN      Summary the result of annotation of domains.
USAGE
	print "$output\n";
}

1;#keep require happy
