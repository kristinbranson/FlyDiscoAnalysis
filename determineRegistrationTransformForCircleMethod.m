function registration_data = determineRegistrationTransformForCircleMethod(expdir, analysis_protocol_folder_path, dataloc_params, registration_params)
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
if isfield(metadata, 'plate') ,
  plate = metadata.plate ;
else
  plate = [] ;
end

% Figure out the actual bowlMarkerType, based on the value in registration
% params, and the plate id, etc
bowlMarkerType = determineBowlMarkerType(registration_params.bowlMarkerType, plate, analysis_protocol_folder_path) ;

% Determine actual maxDistCornerFrac_BowlLabel, possibly by looking it up
% based on the plate ID.
maxDistCornerFrac_BowlLabel = ...
  determineMaxDistCornerFracLabel(registration_params.maxDistCornerFrac_BowlLabel, plate) ;

% Load the bowlMarkerTemplate
if ~endsWith(bowlMarkerType, '.png') ,
  error('Only template-based registration mark types are currently supported.')
end
bowlMarkerTemplate = im2double(imread(bowlMarkerType));  % bowlMarkerType is the name of a .png image, perhaps surprisingly

% Call the core routine to compute the registration transform, and find the
% registration mark(s).
registration_params_cell = filterRegistrationParamsForCircleMethod(registration_params) ;
registration_data = determineRegistrationTransformForCircleMethodCore(registration_params_cell{:}, ...
                                                                      'bkgdImage',bkgdImage, ...
                                                                      'useNormXCorr',true, ...
                                                                      'maxDistCornerFrac_BowlLabel', maxDistCornerFrac_BowlLabel, ...
                                                                      'bowlMarkerTemplate', bowlMarkerTemplate) ;

% Create images illustrating fitting, save to a file as a .png
hfig = [] ;
figpos = [10,10,1600,800] ;
nBowlMarkers = size(registration_data.bowlMarkerPoints,2) ;
saveRegistrationImage(fullfile(expdir, dataloc_params.registrationimagefilestr), ...
                      hfig, ...
                      figpos, ...
                      bkgdImage, ...
                      registration_data.circleCenterX, ...
                      registration_data.circleCenterY, ...
                      registration_data.circleRadius, ...
                      nBowlMarkers, ...
                      registration_data.bowlMarkerPoints, ...
                      registration_params.circleRadius_mm, ...
                      -registration_data.offX, ...
                      -registration_data.offY, ...
                      registration_data.offTheta, ...
                      registration_data.affine, ...
                      registration_params.bowlMarkerPairTheta_true, ...
                      registration_data.scale) ;

