function [success,x,y,featureStrength,filI] = ...
  getNextFeaturePoint(filI, methodcurr, ...
                      nc, nr, ...
                      minFeatureStrengthLow, minFeatureStrengthHigh, ...
                      minTemplateFeatureStrength, ...
                      featureRadius, ...
                      dxGrid, dyGrid)

% find the next strongest feature
[featureStrength,j] = max(filI(:));
[y,x] = ind2sub([nr,nc],j);

% make sure it meets threshold
if strcmpi(methodcurr,'grad2'),
  if featureStrength < minFeatureStrengthHigh,
    success = false ;
    return
  end
else
  if featureStrength < minTemplateFeatureStrength,
    success = false ;
    return
  end
end

% subpixel accuracy
if strcmpi(methodcurr,'grad2'),
  % take box around point and compute weighted average of feature strength
  [box] = padgrab(filI,0,y-featureRadius,y+featureRadius,x-featureRadius,x+featureRadius);
  box = double(box > minFeatureStrengthLow);
  Z = sum(box(:));
  dx = sum(box(:).*dxGrid(:))/Z;
  dy = sum(box(:).*dyGrid(:))/Z;
  x = x + dx;
  y = y + dy;
end

% Zero out region around feature, declare victory
filI = zeroOutDetection(x,y,filI,nc,nr,featureRadius) ;
success = true ;
