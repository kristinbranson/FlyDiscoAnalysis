function result = setPixelsWithinRadiusToNegInf(im, x, y, radius)
% What it says on the tin.
[nr,nc,~] = size(im) ;
result = im ;
ix = round(x) ;
iy = round(y) ;
i1 = max(1,round(ix-radius));  % ALT 2021-03-05: Added outer round() to eliminate warning, should preserve behavior
i2 = min(nc,round(ix+radius));
j1 = max(1,round(iy-radius));
j2 = min(nr,round(iy+radius));
result(j1:j2,i1:i2) = -inf;
