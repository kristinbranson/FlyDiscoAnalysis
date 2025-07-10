% generate the full colormap
n_colors=256;
x=linspace(0,1,n_colors+1)';
x=x(1:end-1);
cmap=oppo3_of_x(x);

% % show it
% figure('color', 'w');
% colormap(cmap);
% colorbar;
% colorbar_axes_h=findobj(gcf,'Tag','Colorbar');
% colorbar_image_h=findobj(colorbar_axes_h,'Tag','TMW_COLORBAR');
% set(colorbar_image_h,'YData',[-180 +180]);
% set(colorbar_axes_h,'YLim',[-180 +180]);
% set(colorbar_axes_h,'YTick',[-180 -90 0 +90 +180]);

% show the colormap a prettier way
theta=(-pi:pi/100:+pi)';
theta=repmat(theta,[1 2]);
r=repmat([0.8 1],[size(theta,1) 1]);
im_index=round(255*(((theta/pi)+1)/2))+1;
im_rgb=ind2rgb(im_index,cmap);
x=r.*cos(theta);
y=r.*sin(theta);
figure('color', 'w');
polar_grid_super_simple();
hold on;
surf(x,y,zeros(size(x)),...
     im_rgb,...
     'EdgeColor','none');
hold off;
text(0,0,'oppo3_of_x',...
     'interpreter','none',...
     'horizontalalignment','center',...
     'verticalalignment','middle');

% make a colormap with inter-color spacings equal
n_colors=256;
cmap=oppo3(n_colors);

% % show the colormap
% figure('color', 'w');
% colormap(cmap);
% colorbar;
% colorbar_axes_h=findobj(gcf,'Tag','Colorbar');
% colorbar_image_h=findobj(colorbar_axes_h,'Tag','TMW_COLORBAR');
% %set(colorbar_image_h,'YData',[-180 +180]);
% set(colorbar_axes_h,'Limits',[-180 +180]);
% set(colorbar_axes_h,'Ticks',[-180 -90 0 +90 +180]);

% show the colormap a prettier way
theta=(-pi:pi/100:+pi)';
theta=repmat(theta,[1 2]);
r=repmat([0.8 1],[size(theta,1) 1]);
im_index=round(255*(((theta/pi)+1)/2))+1;
im_rgb=ind2rgb(im_index,cmap);
x=r.*cos(theta);
y=r.*sin(theta);
figure('color', 'w');
polar_grid_super_simple();
hold on;
surf(x,y,zeros(size(x)),...
     im_rgb,...
     'EdgeColor','none');
hold off;
text(0,0,'oppo3',...
     'interpreter','none',...
     'horizontalalignment','center',...
     'verticalalignment','middle');
