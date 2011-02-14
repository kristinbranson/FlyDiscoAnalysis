#!/usr/bin/perl

`cd /groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis`;
$commit = `git rev-list HEAD --max-count=1`;
print "$commit";
