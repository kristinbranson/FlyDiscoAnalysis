FlyDiscoAnalysis
================

This repo contains code for analyzing FlyDisco experiments.

Entry points that are likely to be useful to end users are `FlyDiscoPipeline()` and
`run_transfero_FlyDiscoPipeline_wrapper_on_experiment_list()`.  See the `help` documentation of those functions for more
details.

To use these functions, start Matlab, `cd` into the root of the repository, and execute the `modpath` command.
This will add the required folders to the Matlab path for the current Matlab session.

This software is designed to be used in conjuction with [Transfero](https://github.com/JaneliaSciComp/transfero), but
using Transfero is in no way a requirement.
