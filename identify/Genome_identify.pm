package Genome_identify;
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use Cwd;

use identify::Identify_Geno_functions;


sub Main
{
	my %options;
	GetOptions(
		"Input:s" => \$options{'Input'},
		"Ref:s" => \$options{'Ref'},
        "Gff:s" => \$options{'Gff'},
		"DB:s" => \$options{'DB'},
		"Outdir:s" => \$options{'Outdir'},
		"Run_type:s" => \$options{'Runtype'},
		"Qsub_queue" => \$options{'Queue'},
		"Retain:s" => \$options{'Retain'},
		"T_Eval:i" => \$options{'T_Eval'},
		"P_Eval:i" => \$options{'P_Eval'},
        "Prefix:s" => \$options{'Prefix'},
		"Process:i" => \$options{'Process'},
        "AA_len:i" => \$options{'AA_len'},
        "Identity:i" => \$options{'Identity'},
        "Max_gap:i" => \$options{'Max_gap'},
	    "Overlap:f" => \$options{'Overlap'},
        "Fa_len:i" => \$options{'Fa_len'},
        "Max_N:f" => \$options{'Max_N'},
    );
	unless($options{'Input'} and $options{'Ref'} and $options{'Gff'})
	{
		return &Usage();
	}
	else
	{
		&Check_opt(\%options);
#print (join "\n",@INC);
        Identify_Geno_functions::Geno_Workflow(\%options);
	}

}

sub Usage
{
	my $usage=<<Usage;

                Usage:  GENOME  -Input <*.pep> -Ref <genome.fa> -DB <validating database> -Gff <*.gff> -Outdir <dir>

                        -Input  Identified genes from proteins or just a batch of genes.
  
                        -Ref    Genome sequences
            
                        -Gff    The gff of input genes.[optional](You can input the gff of whole genome, then you need ensure the GENE_ID is the same as your input)

                        -DB     Validating database[optional]

                        -Outdir    Output directionary

                Options:
                        -Prefix     <str>   the prefix of outputs.[Plant]
                        -T_Eval     <int>   the threshold of E-value of Tblastn.For example, 10 means [e-10].
                        -P_Eval     <int>   the threshold of E-value of blastP. same as T_Eval.[e-10].
                        -Run_type   <Q/P>   how to run those jobs, Q and P represent "qsub" and "process". [P]
                        -Process    <int>   parallel process.[4]
                        -Qsub_queue <str>   if set -Run_type as Q, you can specify the queue to run those jobs. []
                        -Retain     <T/F>   Retain the intermediate result. [F]
                Filter:
                        -AA_len     <int>   For tblastn, minmum AA length, AA_len = Identity * Aligned Length. [disable]
                        -Identity   <int>   For tblastn, minmum identity percent.value range [0,100). [disable]
                        -Max_gap    <int>   For tblastn, Max gap number. [disable]
                        -Overlap    <float> For Checking with gff, value range (0,1], value 1 will be disable the filtration.[0.5]
                        -Fa_len     <int>   For outputing fasta, filtering by the nucl/pep length. [34]
                        -Max_N      <float> For outputing fasta, filtering by N ratio. [0.1]
Usage
	print "$usage\n";
}

sub Check_opt
{
	my $opt = shift;
	my $cwd = getcwd();
	$$opt{'Input'} = abs_path($$opt{'Input'});
    $$opt{'Outdir'}=~s/\/$// if($$opt{'Outdir'});
    $$opt{'Outdir'} = ($$opt{'Outdir'}) ? abs_path("$$opt{'Outdir'}"):"$cwd/Genome_identify";
    $$opt{'Runtype'} ||= "P";
    $$opt{'Retain'} ||= "F";
    $$opt{'Queue'} ||= " ";
    $$opt{'Process'} ||= 4;
    $$opt{'T_Eval'} = &Default_eval($$opt{'T_Eval'});
    $$opt{'P_Eval'} = &Default_eval($$opt{'P_Eval'});
    $$opt{'Fa_len'} ||= 34;
    $$opt{'Max_N'} ||= 0.1;
    $$opt{'Overlap'} ||= 0.5;
}

sub Default_eval
{
    my $eval = shift;
    if($eval)
    {
        return "1e-$eval";
    }
    else
    {
        return "1e-10";
    }
}

1;
#keep require happy
