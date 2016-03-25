package Identify_Pro_functions;
use strict;
use warnings;
use File::Basename;

use identify::Run_jobs;
use identify::Configure;
use identify::Fasta_deal; 
use identify::Protein_identify;

my %Bin;
my ($p_num,$runtype,$queue,$retain);

sub Single_hmm_workflow
{
	my $opt = shift;    #### $opt is a printer;
    Configure::Save_opt($opt);
    ($p_num,$runtype,$queue,$retain) = Configure::Load_opt('Process','Runtype','Queue','Retain');
    &make_dir($opt);
    %Bin = Configure::Bin_Adr();
    my $extract_candidate_fa = &Hmm_iteration($opt);
}

sub make_dir
{
	my $opt = shift;
	unless(-e $$opt{'Outdir'})
	{
		mkdir $$opt{'Outdir'};
		mkdir "$$opt{'Outdir'}/hmmsearch";
		mkdir "$$opt{'Outdir'}/blastp";
	}
    else
    {
        system("rm -r $$opt{'Outdir'}/hmmsearch/*");
        system("rm -r $$opt{'Outdir'}/blastp/*");
    }
}

sub Hmm_iteration
{
    my $opt = shift;
    my ($input,$out_dir) = ($$opt{'Input'},"$$opt{'Outdir'}/hmmsearch");
    my $hmm_file = $$opt{'Hmm_file'}; #### will be replaced
    my $iter_count = 1;
    my $result;
    while($iter_count <= $$opt{'Iter'})
    {
        print "Do interation count : $iter_count\n" if($retain);
        Run_jobs::Print_time();
        my $hmm_result = &Hmmsearch_protein($input,$hmm_file,$iter_count,$out_dir);
        my $hmm_fa = &Extract_candidate_fa($hmm_result,$input);
        Run_jobs::Print_time();
        if($$opt{'DB'})
        {
            my $blastp_result = &Candidate_vs_DB($hmm_fa,$opt);
            my $blastp_fa = &Extract_final_output($blastp_result,$opt);
            Run_jobs::Print_time();
            $hmm_file = &Refresh_hmm($blastp_fa,$iter_count) unless ($iter_count == $$opt{'Iter'});
            Run_jobs::Print_time();
        }
        else
        {
            $hmm_file = &Refresh_hmm($hmm_fa,$iter_count);
        }
        $result = $hmm_fa if($iter_count == $$opt{'Iter'});
        $iter_count ++;
    }
    return $result;
}

sub Hmmsearch_protein
{
    my ($input,$hmm,$count,$out_dir) =@_;
	my $pro_basename = basename $input;
	my $hmm_name = `head -2 $hmm | tail -1`;
	$hmm_name = $1 if($hmm_name=~/NAME\s+(\w+)/);
	$out_dir = $out_dir."/hmm_$hmm_name";
    mkdir $out_dir if(!-e $out_dir);
    my @subfiles;
    if(!-e "$out_dir/$pro_basename.cut/")
    {
        Fasta_deal::Split_fasta_average_size($p_num,$input,$out_dir);
        @subfiles = glob("$out_dir/$pro_basename.cut/*.*");    
    }
    else
    {
        @subfiles = glob("$out_dir/$pro_basename.cut/*.*");
    }
	my $run_shell = "$out_dir/$pro_basename.hmm.sh";
	open OUT,">$run_shell" or die "$!";
	foreach(@subfiles)
	{
		print OUT "$Bin{'hmmsearch'} --noali -E 1e-10 -o $_.log --tblout $_.hmm  $hmm  $_\n";
	}
	close OUT;
	Run_jobs::Run_jobs($run_shell,$runtype,$p_num,$queue);
	system("cat $out_dir/$pro_basename.cut/*.hmm > $out_dir/$pro_basename.$count.hmm");
	return "$out_dir/$pro_basename.$count.hmm";
}

sub Extract_candidate_fa
{
	my ($hmm,$input) = @_;
	my $result_fa = $hmm . ".fa";
	open HMM,"<$hmm" or die "$!";
	my (%id);
	while(<HMM>)
	{
		next if(/#/);
		my @word = split /\s+/;
		$id{$word[0]} = 1;
#print "$word[0]\n";
	}
	close HMM;
    Fasta_deal::Extract_fa_byID(\%id,$input,$result_fa);
	return $result_fa;
}

sub Refresh_hmm
{
    my ($fa,$count) = @_;
    my $msa = "$fa\_$count.msa";
    my $new_hmm = "$fa\_$count.hmm";
    system("$Bin{'clustalw2'} -INFILE=$fa -OUTFILE=$msa -ALIGN -QUIET");
    system("$Bin{'hmmbuild'} -n Iter_$count $new_hmm $msa");
    return $new_hmm;
}

sub Candidate_vs_DB
{
	my ($in_fa,$opt) = @_;
	my $out_dir = $$opt{'Outdir'} . "/blastp";
#system("$Bin{'formatdb'} -i $$opt{'DB'} -p T");
    Fasta_deal::Split_fasta_average_size($p_num,$in_fa,$out_dir);
	my $basename = basename $in_fa;
    my @DB = ($$opt{'DB'});
	my @subfile = glob("$out_dir/$basename.cut/*.*");
	my $run_shell = "$out_dir/$basename.sh";
	my %par;
    %{$par{'blastp'}} = (
            '-evalue' => $$opt{'E_val'},
            '-max_target_seqs' => 1,
            '-outfmt' => 6,
            );
    my $split_dir = $out_dir . "/$basename.cut";
    my @outsub = Run_jobs::Create_blastall(\@subfile,\@DB,\%par,$run_shell,$split_dir);
    Run_jobs::Run_jobs($run_shell,$runtype,$p_num,$queue);
    my $line = join " ",@outsub;
	system("cat $line > $out_dir/$basename.blastp");
	my $outfile = "$out_dir/$basename.blastp";
    return $outfile;
}

sub Extract_final_output
{
	my ($in,$opt) = @_;
	my $out = "$in.fa";
	my %verify;
	open BLAST,"<$in" or die "$!";
	while(<BLAST>)
	{
#next if(/#/);
		chomp;
		my @word = split /\s+/;
		if($word[1]=~/^NBS-LRR/)
		{
			$verify{$word[0]} = 1;
		}
	}
	close BLAST;
    Fasta_deal::Extract_fa_byID(\%verify,$$opt{'Input'},$out);
	return $out;
}

1;
