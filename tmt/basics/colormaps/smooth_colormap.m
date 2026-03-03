function result = smooth_colormap(n_colors, cmap_fn)

% generate a non-smooth colormap on a very fine grid
% calculate the path length, and from that get the phase for each
% point on the fine grid
n_samples=4096;
x=linspace(0,1,n_samples)';
clr=cmap_fn(x);
clr_lab=srgb2lab(clr);
ds=dist_lab(clr_lab(1:end-1,:),clr_lab(2:end,:));
s=[0 ; cumsum(ds)];
phase=s/s(end);  % normalized path length

% make a colormap with inter-color spacings equal
phase_samples=linspace(0,1,n_colors)';
x_samples=interp1(phase,x,phase_samples,'linear');
result=cmap_fn(x_samples);
