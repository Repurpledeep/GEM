package Protein_identify;
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use Cwd;
use Exporter;
our @ISA= qw(Exporter);
our @EXPORT = qw(Main);
our @EXPORT_OK = qw(%options);

use identify::Identify_Pro_functions;

our %options;

sub Main
{
	GetOptions(
		"Input:s" => \$options{'Input'},
		"Hmm_file:s" => \$options{'Hmm_file'},
		"Hmm_list:s" => \$options{'Hmm_list'},
		"DB:s" => \$options{'DB'},
		"Outdir:s" => \$options{'Outdir'},
        "Iter:s" => \$options{'Iter'},
		"Run_type:s" => \$options{'Runtype'},
		"Queue:s" => \$options{'Queue'},
		"Retain:s" => \$options{'Retain'},
		"E_val:i" => \$options{'E_val'},
		"Process:i" => \$options{'Process'},
	);
#	print "$options{'Input'}\n"; #The GetOptions will extract opt from @_
#	exit 1;
    unless($options{'Input'} and ($options{'Hmm_file'} or $options{'Hmm_list'}))
	{
		return &Usage();
	}
	elsif($options{'Hmm_file'} and $options{'Hmm_list'})
	{
		print "\n***Error! only take one of hmm file options***\n";
		return &Usage();
	}
	if($options{'Hmm_file'})
	{
		&Check_and_Set_default(\%options);
		Identify_Pro_functions::Single_hmm_workflow(\%options);
	}
	else
	{
		&Check_and_Set_default(\%options);
		Identify_Pro_functions::Multiple_hmm_workflow(\%options);
	}
}

sub Usage
{
    my $usage =<<USAGE;
        
        Usage:  PROTEIN -Input <*.pep> -Hmm_file <*.hmm> -DB <*.pep> -Outdir
                    
                -Input  Candidate proteins, coding in aa.

                -Hmm_file   A file from hmmsearch for identifying protein with a specific function

                -Hmm_list   A batch of hmm files for identification.[optional]

                -DB Database for refining the result.[optional]

                -Outdir    Output directory. [./Pro_identify]

        Options:
                -Iter   <int>   the times of iteration search with hmmsearch.[3] 
                -E_val  <int>   the threshold of e value.For example, 10 means [e-10].
                -Run_type <Q/P> how to run those jobs, Q and P represent "qsub" and "process". [P]
                -Process <int> parallel process.[8]
                -Qsub_queue if set -Run_type as Q, you may need specify the queue to run those jobs. []
                -Retain <T/F> Retain the intermediate result. [F] 
USAGE
    print "$usage\n";
}

sub Check_and_Set_default
{
	my $opt = shift;
	my $cwd = getcwd();
	$$opt{'Input'} = abs_path($$opt{'Input'});
    $$opt{'DB'} = abs_path($$opt{'DB'});
	$$opt{'Outdir'} = ($$opt{'Outdir'})?abs_path("./$$opt{'Outdir'}"):"$cwd/Pro_identify";
	$$opt{'Runtype'} ||= "P";
	$$opt{'Retain'} ||= "F";
	$$opt{'Queue'} ||= " ";
	$$opt{'Process'} ||= 1;
    $$opt{'Iter'} ||= 5;
	if($$opt{'E_val'})
	{
		$$opt{'E_val'} = "1e-$$opt{'E_val'}";
	}
	else
	{
		$$opt{'E_val'} = 1e-10;
	}
}
1;
#keep require happy
