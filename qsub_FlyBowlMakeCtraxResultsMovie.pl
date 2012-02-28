#!/usr/bin/perl

use strict;
#my $MAXNJOBS = 12;

my $SCRIPTNAME = "FlyBowlMakeCtraxResultsMovie";

my $nargs = $#ARGV + 1;
if($nargs < 1){
    print "Usage: fork_$SCRIPTNAME.pl <expdirlist.txt>\n";
    exit(1);
}

my $ANALYSIS_PROTOCOL = "20120210";

my $SCRIPT = "/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings/$ANALYSIS_PROTOCOL/run_$SCRIPTNAME.sh";

my $PARAMS = "analysis_protocol $ANALYSIS_PROTOCOL";

my $temporary_dir = "/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/temp_$SCRIPTNAME";
`mkdir -p $temporary_dir`;

my $MCR = "/groups/branson/bransonlab/projects/olympiad/MCR/v716";

my $MCR_CACHE_ROOT = "/tmp/mcr_cache_root";

# read in expdirs
my $expdirfile = $ARGV[0];
open(FILE,$expdirfile) or die("Could not open file $expdirfile for reading");

my $njobs = 0;

# loop through each experiment directory
while(my $expdir = <FILE>){
    chomp $expdir;
    
    if(!$expdir){
	next;
    }
    if($expdir =~ /^\s*#/){
	next;
    }

    if(! -d $expdir){
	print "Directory $expdir does not exist\n";
	next;
    }

    $expdir =~ /^(.*)\/([^\/]+)$/;
    my $rootdir = $1;
    my $basename = $2;

    print "*** $basename\n";
    
    # make a name for this job
    my $sgeid = "kb_" . $SCRIPTNAME . "_$basename";
    $sgeid =~ s/\//_/g;
    $sgeid =~ s/\./_/g;
    $sgeid =~ s/\;/_/g;

    # names for temporary script and log file
    my $shfilename = "$temporary_dir/$sgeid" . ".sh";
    my $outfilename = "$temporary_dir/$sgeid" . ".log";

    # create temporary script to be submitted
    write_qsub_sh($shfilename,$expdir,$sgeid);

    # submit command
    my $cmd = qq~qsub -N $sgeid -j y -b y -o $outfilename -cwd $shfilename~;
    print "submitting to cluster: $cmd\n";
    exit;
    system($cmd);
    #system($shfilename);
    $njobs++;

}

print "$njobs jobs started.";
close(FILE);

sub write_qsub_sh {
	my ($shfilename,$expdir,$jobid) = @_;
	
	open(SHFILE,">$shfilename") || die "Cannot write $shfilename";

	print SHFILE qq~#!/bin/bash
# $SCRIPTNAME test script.
# this script will be qsubed
export MCR_CACHE_ROOT=$MCR_CACHE_ROOT.$jobid

echo \$MCR_CACHE_ROOT

$SCRIPT $MCR $expdir $PARAMS

~;
	
#	print SHFILE qq~#delete itself
#rm -f \$0
#~;	
	
	close(SHFILE);
	
	chmod(0755, $shfilename);
}
