function f = figure_of_angular_colormap(cmap, title_string)

theta=(-pi:pi/100:+pi)';
theta=repmat(theta,[1 2]);
r=repmat([0.8 1],[size(theta,1) 1]);
im_index=round(255*(((theta/pi)+1)/2))+1;
im_rgb=ind2rgb(im_index,cmap);
x=r.*cos(theta);
y=r.*sin(theta);
f = figure('color', 'w');
polar_grid_super_simple();
hold on;
surf(x,y,zeros(size(x)),...
     im_rgb,...
     'EdgeColor','none');
hold off;
text(0,0,title_string,...
     'interpreter','none',...
     'horizontalalignment','center',...
     'verticalalignment','middle');
