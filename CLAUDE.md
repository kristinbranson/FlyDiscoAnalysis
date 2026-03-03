This is a project called FlyDiscoAnalysis.  It is a data-processing
pipeline written mostly in Matlab, but with some Python components.
The pipeline takes an input the path of an "experiment directory" or
"expdir", reads raw experimental data in this expdir, and outputs a
bunch of data files that are computed from the raw data.  Each expdir
contains, at a mimimum, a video file named movie.ufmf, and a metadata
file named metadata.xml (with variable capitalizion of the file name).
The movie.ufmf file contains a movie of a bunch of fruit flies in a
chamber, shot from above the chamber.

The pipeline is organized into distinct stages that are run in
sequence.  One of the first stages in the tracking stage, which
extracts the positions and orientations of each fly in each frame.


## FlyTracker calibration

The `n_flies` field inside the parent calibration .mat file
(`flytracker-parent-calibration.mat`) is **not used** during tracking.
Instead, `n_flies` is read from `flytracker-options.txt` in the
analysis protocol folder, and that value overrides whatever is in the
.mat file.  See `core_tracker.m` lines 234-235:

```matlab
if isfield(working_options, 'n_flies') && ~isempty(working_options.n_flies) ,
    pre_fitting_calibration.n_flies = working_options.n_flies ;
```

The same applies to `n_flies_is_max` and `arena_r_mm`.

