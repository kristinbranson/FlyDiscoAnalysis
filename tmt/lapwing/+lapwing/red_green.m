function color_map=red_green(n,q_min,q_max)

% Produces a colormap that maps low values to green, middle values to black,
% and high values to red.  Good for visualizing data that is approximately
% zero-mean, with positive and negative values.  Negative values are red,
% positive are green, and the saturation is a measure of absolute value.
% this version insures that 0 gets mapped to black, q_min to red, and q_max 
% to green

% process args
if nargin<1
  current=colormap;
  n=size(current,1);
end
if nargin==2
  error('Must have 1 or 3 arguments, not 2');
end
if nargin<3
  q_min=-1;
  q_max=+1;
end

% error checking
if q_min>0 || q_max<0
  error('q_min must be nonpositive and q_max must be nonnegative');
end

% make the colormap
ind=(0:(n-1))';
black_index=n*(0-q_min)/(q_max-q_min);
if black_index~=n
  proto=(ind-black_index)/(n-black_index);
else
  proto=zeros(n,1);
end
if black_index~=0
  proto_neg=ind/black_index-1;
  proto(ind<black_index)=proto_neg(ind<black_index);
end
red=-proto.*(proto<0);
green=proto.*(proto>=0);
blue=zeros(n,1);
color_map=[red green blue];
