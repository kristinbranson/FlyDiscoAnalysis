function A = affineTransformMatricesFromOffsetsAndScale(offX,offY,offTheta,scale)
% Compute the affine transforms matrices for a multiple-chamber setup with
% rotation. offX and offY should be chamber_count x 1.  offTheta and scale
% should be scalars.  The result, A, is a 3 x 3 x chamber_count array, with
% A(:,:,k) being the 3x3 affine transform matrix for the kth chamber.  Note
% that a single affine matrix, Ak, is designed to be used like so: 
%
% [x y 1] * A_k => [x_registered y_registered 1]

chamber_count = numel(offX) ;
A = zeros(3,3,chamber_count) ;
Scale = [ scale      0  0 ; ...
              0  scale  0 ; ...
              0      0  1 ] ;
costheta = cos(offTheta) ; 
sintheta = sin(offTheta) ;
Rotation = [  costheta  sintheta  0 ; ...
             -sintheta  costheta  0 ; ...
                     0         0  1 ] ;
for k = 1 : chamber_count ,
  Offset = [       1        0  0 ; ...
                   0        1  0 ; 
             offX(k)  offY(k)  1 ] ;  
  Ak = Offset * Rotation * Scale ;
  A(:,:,k) = Ak ;
end
