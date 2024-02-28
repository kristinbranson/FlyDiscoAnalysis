function [x_reg,y_reg] = registerForSingleChamber(x,y,offX,offY,offTheta,scale)
% [x_reg,y_reg] = register(x,y,offX,offY,offTheta,scale)  
% Given unregistered coordinates <x,y>, add the offset <offX,offY>, rotate
% through angle offTheta, and scale by the *scalar* scale.  <x,y> are assumed to
% be in a Cartesian coordinate system, and the angle offTheta (as is usual)
% specifies a CCW rotation angle.  The resulting registered coordinates
% <x_reg,y_reg> are likewise in a Cartesian coordinate system.  (Yes, all of
% this is pretty much what you'd expect.
sz = size(x);
if numel(sz) ~= numel(size(y)) || ~all(sz == size(y)),
  error('Size of x and y must match');
end
costheta = cos(offTheta); sintheta = sin(offTheta);
X = [x(:)'+offX;y(:)'+offY];
X = [costheta,-sintheta;sintheta,costheta] * X * scale;
x_reg = reshape(X(1,:),sz);
y_reg = reshape(X(2,:),sz);
end  % function
