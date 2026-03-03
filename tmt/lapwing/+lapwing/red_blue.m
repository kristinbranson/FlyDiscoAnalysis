function color_map=red_blue(n)

% Produces a colormap that maps low values to blue, middle values to black,
% and high values to red.  Good for visualizing data that is approximately
% zero-mean, with positive and negative values.  Negative values are red,
% positive are blue, and the saturation is a measure of absolute value.

if nargin<1
  current=colormap;
  n=size(current,1);
end

proto=linspace(-1,1,n)';
red=-proto.*(proto<0);
green=zeros(n,1);
blue=proto.*(proto>=0);
color_map=[red green blue];
