TMT
===

TMT is a library of Matlab functions, centered around the analysis 
of electrophysiology data.


System Requirements
-------------------

Matlab R2013b or later (64-bit)


Installation
------------

1.  Unzip the .zip file to a place of your choosing.

2.  In Matlab, go to File > Set Path... 

3.  Remove any paths for older versions of TMT.

3.  Still in File > Set Path, click on "Add with Subfolders..." and
    select the tmt-release_(whatever) folder that was created when
    you unzipped the .zip file.


Copyright
---------

All files included in TMT present in version 1.09 are copyright Adam
L. Taylor, 1997-2011.  All files after this point, and additions to
existing files, are copyright Howard Hughes Medical Institute,
2011-infinity.


License
-------

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of HHMI nor the names of its contributors may be
  used to endorse or promote products derived from this software
  without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Author
------

[Adam L. Taylor](http://www.janelia.org/people/research-resources-staff/adam-taylor), taylora@hhmi.org  
[Scientific Computing](http://www.janelia.org/research-resources/computing-resources)  
[Janelia Research Campus](http://www.janelia.org)  
[Howard Hughes Medical Institute](http://www.hhmi.org)

[![Picture](/hhmi_janelia_160px.png)](http://www.janelia.org)


Version History
---------------

1.01

Added hds/ directory, the "Hessian direction set" optimization routine.  
Note that it now has a switch so that you can do fixed direction set, 
Powell's method, or my "Hessian direction set" method.

Added help for the function event_triggered().

Deleted event_triggered_double(), since one can just call 
event_triggered() twice.

Deleted fig_triggered_average(), redundant.

Renamed four functions with names like fig_spike_triggered_whatever() to 
fig_event_triggered_whatever, for consistency.

Changed hist_from_data() to more recent version.

Deleted histc_from_data(), so hist_from_data() and histc_from_data() 
weren't out-of-register.  Should probably rewrite histc_from_data() to
be parallel to new hist_from_data().

-----

1.02

Removed cov(), because there's a cov() function built into Matlab that 
does the same thing.  How did I miss that?

Added bipolar(), a colormap function.

Changed to using newer version (2.31) of MATLAB SON Library.  Changed my 
functions that do SON stuff to be compatible.  Ted Brookings did the
modification of load_smr().

Changed to making the user download MATLAB SON Library themselves, because 
Malcolm Lidierth asked us to.

-----

1.03

Moved testing code out of the library, so that adding the library to path
doesn't add any 'modpath.m' files to the path.

Deleted analyze_*.m functions, Z_all_sweeps() from 
data_analysis/impedance_fitting.  Those functions are really not of general 
utility.

Added get_home_dir_name() function.

-----

1.04

Added stg_spectra directory under data_analysis.  This includes code for 
making spectrograms and cohereograms in a way that's tailored for use
with STG rhythms.  Some functions also added in other places to support this.
The code for the plots is similar to what we used for Kris' paper, except we
no longer do the gaussian filtering before calculating the spectra.

-----

1.05

Added function for plotting classic phase diagrams.

Added functions for setting ylims on all axes in a figure.

Added functions mean_circular() and std_circular, for calculating the
mean and standard deviation of angles.

Added function set_axes_size().

Changed print_pdf() to use ghostscript, and to work on UNIX.

Changed print_png() to work on UNIX.

Changed print_tiff() to work on UNIX.

Fixed bug in save_tabular_data() -- didn't open file in text mode.

-----

1.06

Changed some of the spectrogram and coherogram plotting functions to
be more rational.

-----

1.07

Added get_figure_size() function.

Changed times_to_rate_gaussian() so it no longer uses convn().
  Assuming there are not too many spikes, this should make it
  faster most of the time.

Changed print_png(), print_eps(), print_pdf(), print_tiff() to not use
  deprecated -adobecset switch to the built-in print() function.

Deleted show(), because I decided it was too special-purpose.

Deleted chart(), because I decided it was too special-purpose.

Deleted browse(), because I decided it was too special-purpose.

Deleted vcr(), because I decided it was too special-purpose.

Deleted stg_spectra(), because I decided it was too special-purpose.

Deleted impedance_fitting(), because I decided it was too special-purpose.

Deleted load_abf(), replaced it with Michael Sorenson's excellent
  readabf(), because readabf() is better.

Added some documentation to line_wrap().

Added patch_eb(), a simple function to make plotting error bars as patches 
  easier.

Added patch_eb_wrap(), a function for plotting error bars on angular 
  quantities as patches.

Added break_at_wrap_points(), a function that factors out common functionality
  in line_wrap() and patch_eb_wrap().

-----

1.08  (March 7, 2011)

Got rid of readabf(), because had some problems with it.  Put
  load_abf() back in.

Changed all the multitaper code to versions that analyze one signal,
  or a pair of signals, at a time.

Added utility scripts ll and la.

Added get_axes_size()

Added set_axes_size_fixed_center()

Added rate_from_times_and_timeline_gaussian()

Added rate_from_times_and_timeline_simple()

Got rid of plot_wrap()

Got rid of break_at_wrap_points_with_eb(), wasn't used.

Added load_tiff_series()

Added load_multi_image_tiff_file()

-----

1.09 (May 12, 2011)

Fixed bug in coh_mt that lead to NaN's in Cxy_mag getting converted to 1's.

-----

1.10 (August 29, 2011)

Added blue_to_yellow() colormap function.

Added set_axes_position().

Renamed figure_cohereogram() to figure_cohgram().

Renamed figure_spectrogram_amp() to figure_powgram().

Renamed plot_cohereogram() to plot_cohgram().

Renamed plot_spectrogram_amp() to plot_powgram().

Renamed fig_coh_polar() to figure_coh_polar().

Changed figure_cohgram(), plot_cohgram() so that they accept the
  spectrogram in the same units as returned by powgram_mt(), and takes
  an additional argument that specifies the units for plotting.

Changed load_abf() so it now returns the channel names and channel
  units for each channel read in.

Added a function called load_abf_by_channel_name(), which does what
  you would imagine.

Changed patch_eb() so that it divides up the data to make no single
  patch too big.

Fixed patch_eb_wrap() so that it actually returns an array of patch
  object handles.

Added line_eb(), line_eb_wrap() as drop-in replacements for
  patch_eb(), and patch_eb_wrap(), for cases where patches are too
  slow.

Added pow_mt_windowed(), coh_mt_windowed(), that make it more
  convenient to average over windows when calculating spectral
  estimates.

Changed load_named_event_channel_from_smr() to have better error
  messages.

Changed print_pdf() so that it uses the executable gswin64c on 64-bit
  windows machines.

Added tf_mt(), estimates transfer functions of two signals.

Changed coh_mt() to make internal variable names agree with tf_mt(),
  but no user-visible changes.

Added tf_mt_windowed(), analogous to coh_mt_windowed().

Added tfgram_mt(), analogous to cohgram_mt().

Added figure_tf(), figure_coh().

Added a number of functions for reading and writing .tcs files, a file
  format I made up for storing trace data that is easy to import into
  Spike2.

-----

1.11 (December 16, 2011)

Minor changes to several of the spectral plotting functions.

-----

1.12 (March 6, 2012)

Added functions for doing motion-correction of video.

Added functions for doing 'anonymous' saving and loading of .mat files.

-----

1.13 (April 24, 2012)

Minor improvements.

-----

1.14 (November 16, 2012)

Minor improvements.

-----

1.15 (July 28, 2014)

Minor improvements.

-----

1.16 (May 3, 2015)

Fixed bugs with under-specification of character sets in load_abf()
and .tcs reading and writing files.

-----

1.17 (June 29, 2015)

Updated README, added license notice.

