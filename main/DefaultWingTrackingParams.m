function params = DefaultWingTrackingParams()

params = struct;

% background subtraction parameters

% thresholds for wings
params.mindwing_high = 35;
params.mindwing_low = 25;
% thresholds for body
params.mindbody = 100;
% morphology for body
params.radius_dilate_body = 2;
% morphology for wings
params.radius_open_wing = 2;
params.min_wingcc_area = 9;

% maximum distaqnce from a body to wing cc
params.max_wingcc_dist = .25;
% maximum angle to any pixel in the wing
params.max_wingpx_angle = 2.356194490192345;
% obsolete
params.max_wing_angle_otherside = 0.174532925199433;
% minimum non-zero wing angle
params.min_nonzero_wing_angle = 0.174532925199433;
% minimum distance in bins between peaks in the wing angle histogram
params.wing_min_peak_dist_bins = 3;
% minimum fraction of pixels that must fall in the peak bin
params.wing_min_peak_threshold_frac = 0;
% number of bins to histogram the wing angles into
params.nbins_dthetawing = 50;
% minimum fraction of pixels that must fall in the second peak bin, 
% in units of the average fraction per bin
params.wing_peak_min_frac_factor = 2;
% smoothing filter for wing fraction
params.wing_frac_filter = [.25,.5,.25];
% method for fitting wings
params.wing_fit_method = 'peaks';
% obsolete
params.wing_min_prior = .25;
% minimum number of pixels in wings to even try to look for wings
params.min_single_wing_area = 40;
% number of bins to consider for sub-bin accuracy in fitting the wings
params.wing_radius_quadfit_bins = 1;


