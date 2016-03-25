package Run_jobs;
use strict;
use warnings;
use File::Basename;
use POSIX ":sys_wait_h";   ## waitpid;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(Create_blastall Run_jobs);
our @EXPORT_OK = qw(Multi_process);
require identify::Configure;

my %Bin;
sub Create_blastall
{
    %Bin = Configure::Bin_Adr();
    my ($files,$DB,$param,$sh,$out_dir) = @_;
    my @outfile;
    my ($program) = keys %{$param};
#print "$";
    foreach(@{$DB})
    {
#        print SH "####Step_0:Check & make blastdb####\n";
        unless(-e "$_.nhr" || -e "$_.phr")
        {
            if($program eq 'blastx')
            {
                system("$Bin{'formatdb'} -in $_ -dbtype prot\n");
            }
            elsif($program eq 'tblastn')
            {
                system("$Bin{'formatdb'} -in $_ -dbtype nucl\n");
            }
            elsif($program eq 'blastp')
            {
                system("$Bin{'formatdb'} -in $_ -dbtype prot\n");
            }

        }
    }
    open SH,">$sh" or die "$!";
    my $conf = join(" ",%{$param->{$program}});
    for(my $i = 0;$i < @{$DB};$i++)
    {
        my $DB_name = basename $$DB[$i];
        for(my $j = 0;$j < @{$files};$j++)
        {
            my $file_name = basename $$files[$j];
            print SH "$Bin{$program} $conf -db $$DB[$i] -query $$files[$j] -out $out_dir/$file_name\_$DB_name.tabular\n";
            push @outfile,"$out_dir/$file_name\_$DB_name.tabular";
        }
    }
    close SH;
    return @outfile;
}

sub Run_jobs
{
	my ($sh,$type,$p_num,$queue,$retain) = @_;
	if($type eq 'Q')
	{
		system("perl $Bin{'qsub'} $queue --maxjob $p_num $sh");
	}
	elsif($type eq 'P')
	{
        &Multi_process($p_num,$sh,$retain);
    }
	else
	{
		die  "The Run_type is not acceptal to this workflow, Only use\"P\" or \"Q\".\n";
	}
	return 0;
}

sub Merge_file
{
    my ($files,$out_name,$out_dir) = @_;
    open SH,">$out_dir/merge_after_run.sh" or die "$!";
    my $line = join " ",@{$files};
    print SH "cat $line > $out_dir/$out_name.total\n";
    close SH;
    system("sh $out_dir/merge_after_run.sh");
    return "$out_dir/$out_name.total";
}

sub Multi_process
{
    my ($p_num,$sh_file,$retain) = @_;
    $retain ||= "F";
    open LINE,"<$sh_file" or die "$!";
    my @line;
    while(<LINE>)
    {
        chomp;
#        print "$_\n";
        push @line,$_;
    }
    close LINE;
#   print "OK,here $sh_file\n";
    $p_num ||= 1;
#### start multiple process
    my($used_p,$collect_num) = (0 , 0);
    my $collect;
    $SIG{CHLD} = sub { $used_p-- };
    for(my $i = 0 ;$i < @line;$i++)
    {
        my $pid = fork();
        if(!defined($pid))
        {
            die "Error in fork:$!\n";
        }
        if($pid == 0)
        {
            print "Run Paralle jobs:\n$line[$i]\n" if($retain eq "T");
            system ("$line[$i]");
#            sleep(10);
#            print "Finised Paralle jobs:\n$line[$i]\n";
            exit();
        }
        $used_p++;
#        print "use:$used_p\ncount:$i\ncollect:$collect_num\n";
        if( ($i  - $used_p - $collect_num) > 0)
        {
            while( ($collect = waitpid(-1,WNOHANG)) >0)
            {
                $collect_num++;
                print "Done:$collect_num\n" if($retain eq "T");
            }
        }
        do{sleep(1);}until($used_p < $p_num);
    }
    while($collect_num != int(@line))
    {
        while( ($collect = waitpid(-1,WNOHANG)) >0)
        {
            $collect_num++;
            print "Done:$collect_num\n" if($retain eq "T");
        }
    }
}

sub Print_time
{
    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime;
    my $second = time;
    $year += 1900;
    $mon += 1;
    my $datetime = sprintf ("%d-%02d-%02d %02d:%02d:%02d", $year,$mon,$mday,$hour,$min,$sec);
    print $datetime."\n$second\n";
}

1;
