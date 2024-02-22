function bowlMarkerPoints = detectSubimage(im, template, regXY, maxDistCorner_BowlLabel)

% Get image dimensions
[nr,nc,~] = size(im);
r = min(nr,nc);

% compute distance to corners
[xGrid,yGrid] = meshgrid(1:nc,1:nr);
[dxGrid,dyGrid] = meshgrid(-featureRadius:featureRadius,-featureRadius:featureRadius);
distCorner = inf(nr,nc);
corners = [1,1,nc,nc;1,nr,nr,1];
for i = 1:size(corners,2),
  distCorner = min(distCorner, sqrt( (xGrid-corners(1,i)).^2 + (yGrid-corners(2,i)).^2 ));
end

% threshold max distance to some corner
%bowlMarkerTemplate = im2double(imread(bowlMarkerType));  % bowlMarkerType can be a file path...

% Soften template and background
h = fspecial('gaussian',[5,5],1);
template =  imfilter(template,h);
im = imfilter(im,h);

% Scale the template to be on [0,1]
template = template - min(template(:));
template = template / max(template(:));
% If not using normxcorr, want to scale to [-1,+1]
if ~useNormXCorr,
  template = 2*template-1;
end
nRotationsOver180Degrees = 2*nRotations ;
bowlMarkerTemplateFromTheta = cell(1,nRotationsOver180Degrees);
thetas = linspace(0,180,nRotationsOver180Degrees+1);
thetas = thetas(1:end-1);
for i = 1:nRotationsOver180Degrees,
  bowlMarkerTemplateFromTheta{i} = imrotate(template,thetas(i),'bilinear','loose');
end
% compute normalized maximum correlation
filI4_from_theta_index = zeros(nr,nc,nRotationsOver180Degrees) ;
for i = 1:nRotationsOver180Degrees,
  thisBowlMarkerTemplate = bowlMarkerTemplateFromTheta{i} ;
  if useNormXCorr,
    tmp = normxcorr2_padded(thisBowlMarkerTemplate, im, 'replicate') ;
  else
    tmp = imfilter(im,thisBowlMarkerTemplate,            'replicate') ./ ...
      imfilter(im,ones(size(thisBowlMarkerTemplate)),'replicate') ;
  end
  filI4_from_theta_index(:,:,i) = tmp ;
end
filI4 = max(filI4_from_theta_index, [] ,3) ;
methodcurr = 'template';

bowlMarkerIm = filI4; %#ok<NASGU>

% Make non-corners -inf
filI4(distCorner > maxDistCorner_BowlLabel) = -inf;

% make corner with registration mark -inf
if ledindicator
  filI4 = setPixelsWithinRadiusToNegInf(filI4, regXY(1), regXY(2), maxDistCorner_BowlLabel);
end

% find maximum
[success,x,y] = getNextFeaturePoint(filI4, ...
                                    methodcurr, ...
                                    nc, nr, ...
                                    minFeatureStrengthLow, minFeatureStrengthHigh, ...
                                    minTemplateFeatureStrength, ...
                                    featureRadius, ...
                                    dxGrid, dyGrid);
if success,
  bowlMarkerPoints = [x y]' ;
else
  error('Could not detect bowl marker');
end
end  % function








