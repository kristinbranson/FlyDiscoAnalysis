function [success,x,y] = ...
  getFeaturePoint(filI, ...
                  minTemplateFeatureStrength)

% find the strongest feature
[nc,nr,~] = size(filI) ;
[featureStrength,j] = max(filI(:));
[y,x] = ind2sub([nr,nc],j);

% make sure it meets threshold
success = (featureStrength >= minTemplateFeatureStrength) ;