#!/usr/bin/perl

use strict;
use Cwd;

my $ANALYSISDIR = "/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis";
my $CTRAXDIR = "/groups/branson/bransonlab/projects/Ctrax/Ctrax";
my $CTRAXSETTINGSDIR = "/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax";
my $OUTFILESTR = "analysis_protocol.txt";


my $nargs = $#ARGV + 1;
if($nargs < 1){
    print "Usage: analysis_protocol.pl <expdir>\n";
    exit(1);
}
my $expdir = $ARGV[0];
if(! -d $expdir){
    print "Directory $expdir does not exist\n";
    exit(1);
}

# get git hash tag
my $dir0 = getcwd;
chdir $ANALYSISDIR;
my $analysis_version = `git rev-list HEAD | head -n 1`;
chomp $analysis_version;
chdir $dir0;
print "analysis version: $analysis_version\n";

# get settings parameters
my $line = `ls -l $ANALYSISDIR/settings/current`;
$line =~ /-> (.*\/)?([^\/]*)\/?$/;
my $analysis_settings = $2;
chomp $analysis_settings;
print "analysis settings: $analysis_settings\n";

my $ctrax_version = "";

# get Ctrax version
open(FILE,"$CTRAXDIR/version.py");
while($line = <FILE>){
    chomp $line;
    $line =~ /^__version__ = "([^"]+)"$/;
    if($1){
	$ctrax_version = $1;
	last;
    }
}
close(FILE);
print "ctrax version: $ctrax_version\n";

# get Ctrax settings parameters
my $line = `ls -l $CTRAXSETTINGSDIR/current`;
$line =~ /-> (.*\/)?([^\/]*)\/?$/;
my $ctrax_settings = $2;
chomp $ctrax_settings;
print "ctrax settings: $ctrax_settings\n";

my $analysis_protocol = "ctrax_version:$ctrax_version,ctrax_settings:$ctrax_settings:$ctrax_settings,analysis_settings:$analysis_settings,analysis_version:$analysis_version";

print "analysis_protocol:\n$analysis_protocol\n";
open(FILE,">$expdir/$OUTFILESTR");
print FILE "$analysis_protocol\n";
close(FILE);
