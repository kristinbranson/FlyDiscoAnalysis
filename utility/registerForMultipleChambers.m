function [x_reg,y_reg,closest_chamber_index] = registerForMultipleChambers(x,y,offX,offY,offTheta,scale)
% [x_reg,y_reg] = registerForMultipleChambers(x, y, offX, offY, offTheta, scale) 
% 
% Given unregistered coordinates <x,y>, add the offset <offX,offY>, translate
% by <offX,offY>, rotate by the *scalar* offTheta, and scale by the *scalar*
% scale.  <x,y> are assumed to be in a Cartesian coordinate system.  x and y
% should be of the same shape. offX and offY should be chamber_count x 1
% arrays, with the offset for the jth chamber in the jth row.  Note this
% function does not support rotation.

% Throw a more-scrutable error if x and y differ in shape
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

% Coerce offX and offY to col vectors
offX_col = offX(:) ;
offY_col = offY(:) ;

% For each <x,y>, determine which chamber it belongs to
distance_to_center = hypot(x_row+offX_col, y_row+offY_col) ;  % chamber_count x n
[~, closest_chamber_index_as_row] = min(distance_to_center, [], 1) ;

% Get the center of the nearest chamber to each point
n = numel(closest_chamber_index_as_row) ;
offXFull = reshape(offX_col(closest_chamber_index_as_row), [1 n]) ;  % 1 x n, x-offset of nearest chamber center
offYFull = reshape(offY_col(closest_chamber_index_as_row), [1 n]) ;  % 1 x n, y-offset of nearest chamber center

% Compute the rotation matrix
% Note that rotation is, in principle, about the mean of all the per-chamber origins, but this still
% works.
costheta = cos(offTheta); sintheta = sin(offTheta);
R = [costheta,-sintheta;sintheta,costheta] ;

% Do the transform
Xt = [ x_row+offXFull ; ...
       y_row+offYFull ] ;
Y = R * Xt * scale ;
x_reg_row = Y(1,:) ;
y_reg_row = Y(2,:) ;

% Coerce the outputs to be the same shape as the original x and y
if were_x_and_y_both_rows ,
  % do nothing
  x_reg = x_reg_row ;
  y_reg = y_reg_row ;
  closest_chamber_index = closest_chamber_index_as_row ;
else
  x_reg = reshape(x_reg_row, sz) ;
  y_reg = reshape(y_reg_row, sz) ;
  closest_chamber_index = reshape(closest_chamber_index_as_row, sz) ;
end
