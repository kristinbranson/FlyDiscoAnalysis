function [x_reg,y_reg] = registerForMultipleChambers(x,y,offX,offY,offTheta,scale)
% [x_reg,y_reg] = registerForMultipleChambers(x, y, offX, offY, offTheta, scale) 
% 
% Given unregistered coordinates <x,y>, add the offset <offX,offY>, translate
% by <offX,offY>, rotate by the *scalar* offTheta, and scale by the *scalar*
% scale.  <x,y> are assumed to be in a Cartesian coordinate system.  x and y
% should be of the same shape. offX and offY should be chamber_count x 1
% arrays, with the offset for the jth chamber in the jth row.  Note this
% function does not support rotation.

sz = size(x);
if numel(sz) ~= numel(size(y)) || ~all(sz == size(y)),
  error('Size of x and y must match');
end

% Coerce x,y to row vectors
were_x_and_y_both_rows = isrow(x) && isrow(y) ;
if were_x_and_y_both_rows ,
  % do nothing
  x_row = x ;
  y_row = y ;
else
  x_row = x(:)' ;
  y_row = y(:)' ;
end

% For each <x,y>, determine which chamber it belongs to
distance_to_center = hypot(x_row+offX, y_row+offY) ;  % chamber_count x n
[~, closest_chamber_index] = min(distance_to_center, [], 1) ;

% Get the center of the nearest chamber to each point
n = numel(closest_chamber_index) ;
offXFull = reshape(offX(closest_chamber_index), [1 n]) ;  % 1 x n, x-offset of nearest chamber center
offYFull = reshape(offY(closest_chamber_index), [1 n]) ;  % 1 x n, y-offset of nearest chamber center

% Compute the rotation matrix
% Note that rotation is, in principle, about the mean of all the per-chamber origins, but this still
% works.
costheta = cos(offTheta); sintheta = sin(offTheta);
R = [costheta,-sintheta;sintheta,costheta] ;

% Do the transform
Xt = [ x+offXFull ; ...
       y+offYFull ] ;
Y = R * Xt * scale ;
x_reg_row = Y(1,:) ;
y_reg_row = Y(2,:) ;

% Coerce the outputs to be the same shape as the original x and y
if were_x_and_y_both_rows ,
  % do nothing
  x_reg = x_reg_row ;
  y_reg = y_reg_row ;
else
  x_reg = reshape(x_reg_row, sz) ;
  y_reg = rehapse(y_reg_row, sz) ;
end
