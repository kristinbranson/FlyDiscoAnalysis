#!/usr/bin/perl

use strict;
my $MAXNJOBS = 12;

my $nargs = $#ARGV + 1;
if($nargs < 1){
    print "Usage: fork_FlyBowlComputePerFrameFeatures.pl <expdirlist.txt>\n";
    exit(1);
}

my $ANALYSIS_PROTOCOL = "20110407";

my $SCRIPT = "/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings/$ANALYSIS_PROTOCOL/run_FlyBowlComputePerFrameFeatures.sh";

my $PARAMS = "analysis_protocol $ANALYSIS_PROTOCOL";

my $temporary_dir = "/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/temp_compute_stats";
`mkdir -p $temporary_dir`;

my $MCR = "/groups/branson/bransonlab/projects/olympiad/MCR/v714";

my $MCR_CACHE_ROOT = "/tmp/mcr_cache_root";

# read in expdirs
my $expdirfile = $ARGV[0];
open(FILE,$expdirfile) or die("Could not open file $expdirfile for reading");

my $njobs = 0;

# loop through each experiment directory
my @childs = ();
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
    my $sgeid = "kb_flybowlcomputeperframefeatures_$basename";
    $sgeid =~ s/\//_/g;
    $sgeid =~ s/\./_/g;
    $sgeid =~ s/\;/_/g;

    # names for temporary script and log file
    my $shfilename = "$temporary_dir/$sgeid" . ".sh";
    my $outfilename = "$temporary_dir/$sgeid" . ".log";

    # create temporary script to be submitted
    write_qsub_sh($shfilename,$expdir,$sgeid);

    # submit command
    my $cmd = "$shfilename 2>&1 > $outfilename";

    my $pid = fork();
    if($pid){
	$njobs++;
	push(@childs, $pid);
    }
    else{
	print "job $njobs: $cmd\n";
	#sleep(10);
	`$cmd`;
	$njobs--;
	exit(0);
    }


    #system($cmd);
    #system($shfilename);


    if($njobs >= $MAXNJOBS){
	print "njobs = $njobs\n";
	my $child = wait();
	print "child $child finished\n";
	$njobs--;
    }

}

close(FILE);



foreach (@childs) {
    waitpid($_, 0);
    print "child $_ finished\n";
}


sub write_qsub_sh {
	my ($shfilename,$expdir,$jobid) = @_;
	
	open(SHFILE,">$shfilename") || die "Cannot write $shfilename";

	print SHFILE qq~#!/bin/bash
# FlyBowlComputePerFrameFeatures test script.
# this script will be qsubed
export MCR_CACHE_ROOT=$MCR_CACHE_ROOT.$jobid

$SCRIPT $MCR $expdir $PARAMS

~;
	
#	print SHFILE qq~#delete itself
#rm -f \$0
#~;	
	
	close(SHFILE);
	
	chmod(0755, $shfilename);
}
