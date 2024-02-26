function registration_data = determineRegistrationTransformForFlyTrackerMethod(expdir, analysis_protocol_folder_path, dataloc_params, registration_params)
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

% % Call the core routine to compute the registration transform, and find the
% % registration mark(s).
% registration_params_cell = marshallRegistrationParams(registration_params, expdir, dataloc_params.registrationimagefilestr) ;
% registration_data = determineRegistrationTransformForCircleMethodCore(registration_params_cell{:}, ...
%                                                                       'bkgdImage',bkgdImage, ...
%                                                                       'useNormXCorr',true, ...
%                                                                       'bowlMarkerType', bowlMarkerType, ...
%                                                                       'maxDistCornerFrac_BowlLabel', maxDistCornerFrac_BowlLabel) ;
%   % bowlMarkerType and maxDistCornerFrac_BowlLabel will override values in registration_params_cell
useNormXCorr = true ;

% Get basic dimensions
[nr,nc,~] = size(bkgdImage);
r = min(nr,nc);



%
% Determine where the chamber boundary is
%

% Only circle and flytracker methods are supported in FlyDisco
if ~(strcmpi(method, 'circle') || strcmpi(method, 'flytracker')) ,
  error('Unsupported method for chamber boundary detection: %s',method);  
end
% Depending on the circle image type, compute an image (circleim) that will be used for
% circle detection.
if strcmpi(circleImageType,'raw_whiteedge'),
  circleim = bkgdImage >= circleImageThresh;
elseif strcmpi(circleImageType,'raw_blackedge'),
  circleim = bkgdImage <= circleImageThresh;
elseif strcmpi(circleImageType,'grad'),
  gradI = [diff(bkgdImage,1,1).^2;zeros(1,nc)] + [diff(bkgdImage,1,2).^2,zeros(nr,1)];
  circleim = sqrt(gradI) >= circleImageThresh;
elseif strcmpi(circleImageType,'canny'),
  circleim = bkgdImage;
else
  error('Unknown circleImageType %s',circleImageType);
end
% Call the core circle-detection routine
binedgesx = linspace(circleXLim(1)*nc,circleXLim(2)*nc,circleNXTry+1);
bincentersy = linspace(circleYLim(1)*nr,circleYLim(2)*nr,circleNYTry);
bincentersr = linspace(circleRLim(1)*min(nc,nr),circleRLim(2)*min(nc,nr),circleNRTry);
[circleRadius,circleCenterX,circleCenterY,featureStrengths,circleDetectParams] = ...
  detectcircles(circleim,...
                'cannythresh',circleCannyThresh,'cannysigma',circleCannySigma,...
                'binedgesa',binedgesx,'bincentersb',bincentersy,'bincentersr',bincentersr,...
                'maxncircles',1,'doedgedetect',strcmpi(circleImageType,'canny'));



%
% Find the registration mark
%

% compute distance to corners
if nBowlMarkers>1 ,
  error('Only zero or one registration mark is currently supported') ;
end
if nBowlMarkers == 0 ,
  bowlMarkerPoints = zeros(2,0);
else
  if ~endsWith(bowlMarkerType, '.png') ,
    error('Only template-based registration mark types are currently supported.')
  end
  bowlMarkerTemplate = im2double(imread(bowlMarkerType));  % bowlMarkerType is the name of a .png image, perhaps surprisingly
  corner_radius = maxDistCornerFrac_BowlLabel * r ;
  bowlMarkerPoints = ...
    findTemplateMatchWithPossibleRotation(bkgdImage, bowlMarkerTemplate, minTemplateFeatureStrength, nRotations, useNormXCorr, corner_radius) ;
end



%
% Define the paramters of the transform used to register tracks.
%

% The origin is the center of the chamber circle.
originX = circleCenterX;
originY = circleCenterY;

% The offset is the vector you add to unregistered points to translate them
% s.t. points at the origin have registered coordinates of (0,0)/
offX = -originX;
offY = -originY;

% Determine the angle of the registration mark
if nBowlMarkers > 0 ,
  bowlMarkerPoint = mean(bowlMarkerPoints, 2, 'omitnan') ;
  bowlMarkerTheta = atan2(bowlMarkerPoint(2)-originY,bowlMarkerPoint(1)-originX);
else
  % Pretend the registration mark is in the upper-left corner of the image
  bowlMarkerTheta = atan2(1-originY,1-originX);  
end

% For FlyDisco, we assume that no rotation is needed to bring things into
% register.  In the FlyBowl era, there was often a n*90 degree rotation needed
% to bring things into register.
offTheta = 0;

% Compute the scale factor
scale = circleRadius_mm / circleRadius;  % mm/pixel



%
% Package everything up into the registration struct, save that to disk
% 
registerfn = @(x,y)(registerForSingleAffine(x,y,offX,offY,offTheta,scale)) ;
affine = affineTransformMatrixFromOffsetsAndScale(offX,offY,offTheta,scale);
registration_data = ...
  struct('offX',offX,...
         'offY',offY,...
         'offTheta',offTheta,...
         'scale',scale,...
         'bowlMarkerTheta',bowlMarkerTheta,...
         'bkgdImage',bkgdImage,...
         'featureStrengths',featureStrengths,...
         'affine',affine, ...
         'bowlMarkerPoints', bowlMarkerPoints, ...
         'registerfn', registerfn, ...
         'circleCenterX', circleCenterX, ...
         'circleCenterY', circleCenterY, ...
         'circleRadius', circleRadius) ;
% Want registration to be a scalar struct, but circleDetectParams is a cell
% array.  So easier to tack it on here.
registration_data.circleDetectParams = circleDetectParams ;  

% create images illustrating fitting
saveRegistrationImage(imsavename, hfig, figpos, bkgdImage, circleCenterX, circleCenterY, circleRadius, ...
                      nBowlMarkers, bowlMarkerPoints, circleRadius_mm, originX, originY, offTheta, affine, ...
                      bowlMarkerPairTheta_true, scale) ;






