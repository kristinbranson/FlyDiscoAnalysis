function registration_data = determineRegistrationTransformForFlyTrackerMethodCore(varargin)
% The core routine for fitting the chamber border, computing the registration
% transform, and finding the registration mark, when using the 'circle'
% registration method.
%
% This is a pure function, unless debug is true.  Then it makes some plots but
% it otherwise pure.

% parse inputs
[bkgdImage,...
 nRotations,...
 minTemplateFeatureStrength,...
 maxDistCornerFrac_BowlLabel,...
 bowlMarkerTemplateOrEmpty,...
 isInDebugMode,...
 useNormXCorr, ...
 flytrackerCalibration] = ...
  myparse(varargin,...
  'bkgdImage',[],...
  'nRotations',20,...
  'minTemplateFeatureStrength',.92,...
  'maxDistCornerFrac_BowlLabel',.17,...
  'bowlMarkerTemplateOrEmpty',[],...
  'debug',false,...
  'useNormXCorr',false, ...
  'flytrackerCalibration', []);


% Error if required args not supplied
if isempty(bkgdImage) ,
  error('No background image supplied to %s()', mfilename()) ;
end
if isempty(flytrackerCalibration) ,
  error('No FlyTracker output calibration supplied to %s()', mfilename()) ;
end

% Get basic dimensions
[nr,nc,~] = size(bkgdImage);



%
% Determine where the chamber boundary is
%

% Read the centroids, common radius out of the FT calibration file
centroids = flytrackerCalibration.centroids ;  % chamberCount x 2, in pels
r = flytrackerCalibration.r ;  % pels, scalar for circular chambers
if isempty(r) ,
  error('The FlyTracker calibration file "r" field is empty.  Currently, only circular chambers are supported.') ;
end
% chamberCount = size(centroids,1) ;
circleRadius = r ;  % scalar, b/c all arenas have to be the same size
circleCenterX = centroids(:,2) ;  % chamberCount x 1
circleCenterY = centroids(:,1) ;  % chamberCount x 1



%
% Find the registration mark
%

% compute distance to corners
if isempty(bowlMarkerTemplateOrEmpty) ,
  % No bowl marker need be detected.
  bowlMarkerPoints = zeros(2,0);
else
  bowlMarkerTemplate = bowlMarkerTemplateOrEmpty ;
  corner_radius = maxDistCornerFrac_BowlLabel * min(nr,nc) ;
  bowlMarkerPoints = ...
    findTemplateMatchWithPossibleRotation(bkgdImage, bowlMarkerTemplate, minTemplateFeatureStrength, nRotations, useNormXCorr, corner_radius) ;
end



%
% Define the paramters of the transform used to register tracks.
%

% The origin is the center of the chamber circle.
originX = circleCenterX;  % chamberCount x 1, in pels
originY = circleCenterY;  % chamberCount x 1, in pels

% The offset is the vector you add to unregistered points to translate them
% s.t. points at the origin have registered coordinates of (0,0)/
offX = -originX;  % chamberCount x 1, in pels
offY = -originY;  % chamberCount x 1, in pels

% Determine the angle of the registration mark
meanCentroid = mean(centroids,1) ;  % 1 x 2, in pels
if ~isempty(bowlMarkerTemplateOrEmpty) ,
  bowlMarkerPoint = mean(bowlMarkerPoints, 2, 'omitnan') ;
  bowlMarkerTheta = atan2(bowlMarkerPoint(2)-meanCentroid(2),bowlMarkerPoint(1)-meanCentroid(1));
else
  % Pretend the registration mark is in the upper-left corner of the image
  bowlMarkerTheta = atan2(1-meanCentroid(2),1-meanCentroid(1));  
end

% For FlyDisco, we assume that no rotation is needed to bring things into
% register.  In the FlyBowl era, there was often a n*90 degree rotation needed
% to bring things into register.
offTheta = 0;   % radians

% Compute the scale factor
circleRadius_mm = flytrackerCalibration.arena_r_mm ;  % mm
scale = circleRadius_mm / r;  % mm/pixel



%
% Package everything up into the registration struct, save that to disk
% 
registerfn = @(x,y)(registerForMultipleChambers(x,y,offX,offY,offTheta,scale)) ;
affines = affineTransformMatricesFromOffsetsAndScale(offX,offY,offTheta,scale) ;
registration_data = ...
  struct('offX',offX,...
         'offY',offY,...
         'offTheta',offTheta,...
         'scale',scale,...
         'bowlMarkerTheta',bowlMarkerTheta,...
         'bkgdImage',bkgdImage,...
         'affine',affines, ...
         'bowlMarkerPoints', bowlMarkerPoints, ...
         'registerfn', registerfn, ...
         'circleCenterX', circleCenterX, ...
         'circleCenterY', circleCenterY, ...
         'circleRadius', circleRadius) ;

% % Make debug plot(s)
% iscircle = true ;
% makeRegistrationDebugPlots(isInDebugMode, iscircle, bkgdImage, circleRadius_mm, pairDist_mm, originX, originY, offTheta, ...
%                            circleCenterX, circleCenterY, circleRadius, [], bowlMarkerPoints, ...
%                            bowlMarkerPairTheta_true, markerPairAngle_true) ;





