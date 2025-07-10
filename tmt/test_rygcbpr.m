% generate the full colormap
n_colors=256;
x=linspace(0,1,n_colors+1)';
x=x(1:end-1);
cmap0=rygcbpr_of_x(x);

% show the colormap a pretty way
f0 = figure_of_angular_colormap(cmap0, 'rygcbpr_of_x') ;

% make a colormap with inter-color spacings equal
n_colors=256;
cmap=smooth_colormap(n_colors, @rygcbpr_of_x);

% show the colormap a pretty way
f = figure_of_angular_colormap(cmap, 'rygcbpr') ;
