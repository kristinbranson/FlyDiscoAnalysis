function result = findTemplateMatchWithPossibleRotation(im, raw_template, minTemplateFeatureStrength, nRotations, useNormXCorr, corner_radius, excluded_xyrs)
% Find a place in the image im that matches the raw_template.  Result is a 2x1
% col vector containing x,y coords.

% Handle args
if ~exist('corner_radius', 'var') || isempty(corner_radius) ,
  corner_radius = inf ;
end
if ~exist('excluded_xyrs', 'var') || isempty(excluded_xyrs) ,
  excluded_xyrs = zeros(0,3) ;
end

% Get image dimensions
[nr,nc,~] = size(im);

% Soften template and background
kernel = fspecial('gaussian',[5,5],1);
template = imfilter(raw_template,kernel);
softIm = imfilter(im,kernel);

% Scale the template to be on [0,1]
template = template - min(template(:));
template = template / max(template(:));

% If not using normxcorr, want to scale to [-1,+1]
if ~useNormXCorr,
  template = 2*template-1;
end

% Rotate the template, to match any rotation
nRotationsOver180Degrees = 2*nRotations ;
bowlMarkerTemplateFromTheta = cell(1,nRotationsOver180Degrees);
thetas = linspace(0,180,nRotationsOver180Degrees+1);
thetas = thetas(1:end-1);
for i = 1:nRotationsOver180Degrees,
  bowlMarkerTemplateFromTheta{i} = imrotate(template,thetas(i),'bilinear','loose');
end
% compute normalized maximum correlation
match_image_from_theta_index = zeros(nr,nc,nRotationsOver180Degrees) ;
for i = 1:nRotationsOver180Degrees,
  thisBowlMarkerTemplate = bowlMarkerTemplateFromTheta{i} ;
  if useNormXCorr,
    match_image = normxcorr2_padded(thisBowlMarkerTemplate, softIm, 'replicate') ;
  else
    match_image = imfilter(softIm,thisBowlMarkerTemplate,            'replicate') ./ ...
                  imfilter(softIm,ones(size(thisBowlMarkerTemplate)),'replicate') ;
  end
  match_image_from_theta_index(:,:,i) = match_image ;
end
rawTemplateMatchImage = max(match_image_from_theta_index, [] ,3) ;

% compute distance to corners
[xGrid,yGrid] = meshgrid(1:nc,1:nr);
distanceToCornerImage = inf(nr,nc);
corners = [1,1,nc,nc;1,nr,nr,1];
for i = 1:size(corners,2),
  distanceToCornerImage = min(distanceToCornerImage, sqrt( (xGrid-corners(1,i)).^2 + (yGrid-corners(2,i)).^2 ));
end

% Make non-corners -inf
templateMatchImage = rawTemplateMatchImage ;
templateMatchImage(distanceToCornerImage > corner_radius) = -inf;

% Make points within exclustion zones -inf
exclusion_zone_count = size(excluded_xyrs, 1) ;
for zone_index = 1 :  exclusion_zone_count ,
  excluded_xyr = excluded_xyrs(zone_index,:) ;
  templateMatchImage = setPixelsWithinRadiusToNegInf(templateMatchImage, excluded_xyr(1), excluded_xyr(2), excluded_xyr(3)) ;
end

% find maximum
[success,x,y] = getFeaturePoint(templateMatchImage, ...
                                minTemplateFeatureStrength) ;
if success,
  result = [x y]' ;
else
  error('Could not find a good engouh match to template');
end
