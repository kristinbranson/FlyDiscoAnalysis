function registration_data = determineRegistrationTransform(expdir, dataloc_params, registration_params)
% Determine the registration transform.  Along the way, detect the
% registration mark(s), even though they are currently not used for
% registration.

% Load the background file output by FlyTracker
if isfield(dataloc_params,'flytrackerbgstr')
  flytrackerbgfile = fullfile(expdir,dataloc_params.flytrackerbgstr);
  load(flytrackerbgfile,'bg');
  bkgdImage = 255*bg.bg_mean;  % double image, each el on [0,255], not necessarily integers
else
  % Note: We really do want to error if this is missing.
  % We've been bitten by this not being what we thought it was.
  error('dataloc_params is missing field flytrackerbgstr, which is required');
end

% Load metadata
metadata = collect_metadata(expdir, dataloc_params.metadatafilestr) ;
plate = metadata.plate ;

% Figure out the actual bowlMarkerType, based on the value in registration
% params, and the plate id, etc
bowlMarkerType = determineBowlMarkerType(registration_params.bowlMarkerType, plate, analysis_protocol_folder_path) ;

% Determine actual maxDistCornerFrac_BowlLabel, possibly by looking it up
% based on the plate ID.
maxDistCornerFrac_BowlLabel = ...
  determineMaxDistCornerFracLabel(registration_params.maxDistCornerFrac_BowlLabel, plate) ;

% Call the core routine to compute the registration transform, and find the
% registration mark(s).
registration_params_cell = marshallRegistrationParams(registration_params, expdir, dataloc_params.registrationimagefilestr) ;
registration_data = determineRegistrationTransformCore(registration_params_cell{:}, ...
                                                       'bkgdImage',bkgdImage, ...
                                                       'useNormXCorr',true, ...
                                                       'bowlMarkerType', bowlMarkerType, ...
                                                       'maxDistCornerFrac_BowlLabel', maxDistCornerFrac_BowlLabel) ;
  % bowlMarkerType and maxDistCornerFrac_BowlLabel will override values in registration_params_cell
