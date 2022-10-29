FlyDiscoAnalysis
================

This repo contains code for analyzing FlyDisco experiments.  It includes Goldblum,
a service for automatically tranferring experiment folders off of a set of rig computers, and 
the running the FlyDiscoAnalysis pipeline on them.

Entry points that are likely to be useful to end users are `FlyDiscoPipeline()` and
`goldblum_analyze_experiment_folders()`.  See the `help` documentation of those functions for more
details.

To use these functions, start Matlab, `cd` into the root of the repository, and execute the `modpath` command.
This will add the required folders to the Matlab path for the current Matlab session.

The [Goldblum README](goldblum/README.md) provides more details on what the Goldblum service is,
and how to set it up.
