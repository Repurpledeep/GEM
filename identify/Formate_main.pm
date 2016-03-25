package Formate_main;
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';

use identify::Fasta_deal;

sub Main
{
    if( int(@_) < 2)
    {
        &Main_Usage();
    }
    else
    {
        shift @_; ## delete the main "Formate";
        my $sub = shift @_;
        print "$sub\n";
        if($sub eq "GFF2PEP")
        {
            &GFF_PEP_main(@_);
        }
        elsif($sub eq "SPLIT")
        {
        }
        else
        {
            &Main_Usage();
        }
    }
}

sub Main_Usage
{
    my $main_usage =<<USAGE;

        Formate Usage:
                
                GFF2PEP     Input the gff and ref sequences to extract protein sequences.
                SPLIT       Split the fasta file into several sub-files[developing]. 

USAGE
    print $main_usage;
}
 
sub GFF_PEP_main
{
    my %opt;
    GetOptions(
            "Gff:s" => \$opt{'Gff'},
            "Ref:s" => \$opt{'Ref'},
            "Output:s" => \$opt{'Output'},
            "M_len:i" => \$opt{'M_len'},
            "N_ratio:f" => \$opt{'N_ration'},
            );
    my $usage =<<USAGE_1;

        GFF2PEP     -Gff <*.gff> -Ref <*.fa> -Output <*.pep>
                
                -Gff        Gff file.
                -Ref        Reference sequences.
                -Output     Output file, null will output to STD.[optional]
        
        Options:
                
                -M_len      Output threshold of minimum aa length.[34]
                -N_ratio    Output threshold of maximum N ratio,value range(0,1]. [0.1]

USAGE_1
    if(!$opt{'Gff'} || !$opt{'Ref'})
    {
        print "Error:Lack of opt.\n$usage";
        exit 1;
    }
    elsif(! (-e $opt{'Gff'}) || ! (-e $opt{'Ref'}))
    {
        print "Error: No such file.\n$usage";
        exit 1;
    }
    else
    {
        $opt{'Gff'} = abs_path($opt{'Gff'});
        $opt{'Ref'} = abs_path($opt{'Ref'});
        $opt{'M_len'} ||= 34;
        $opt{'N_ratio'} = 0.1 if(!$opt{'N_ratio'} || $opt{'N_ratio'} > 1);
        my %pep = Fasta_deal::Gff_to_Pep($opt{'Gff'},$opt{'Ref'});
        if($opt{'Output'})
        {
            Fasta_deal::Print_fasta(\%pep,$opt{'Output'},$opt{'M_len'},$opt{'N_ratio'});
        }
        else
        {
            Fasta_deal::Print_fasta_STD(\%pep,,$opt{'M_len'},$opt{'N_ratio'});
        }
    }
}

1; ### keep require happy
