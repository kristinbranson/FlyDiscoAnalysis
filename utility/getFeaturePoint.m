function [success,x,y] = ...
  getFeaturePoint(templateMatchImage, ...
                  minTemplateMatchStrength)

% find the strongest feature
assert(ismatrix(templateMatchImage), 'Internal error: templateMatchImage in getFeaturePoint() must be a matrix') ;
[nr, nc] = size(templateMatchImage) ;
[featureStrength, j] = max(templateMatchImage(:)) ;
[y, x] = ind2sub([nr nc], j) ;

% make sure it meets threshold
success = (featureStrength >= minTemplateMatchStrength) ;

end  % function
