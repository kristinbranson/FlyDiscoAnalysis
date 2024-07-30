function saveLedDetectionImage(imsavename, bkgdImage, ledXY, templateShapeXY, doesYAxisPointUp)

if isempty(imsavename) 
  return
end

figpos = [10,10,800,800];
hfig = figure();
set(hfig,'Units','pixels','Position',figpos);
nimsplot = 1;
hax = createsubplots(1,nimsplot,0.05,hfig);
 
% plot background image in jet colormap
imagesc(bkgdImage,'parent',hax(1),[0,255]);
axis(hax(1),'image');
if doesYAxisPointUp ,
  axis(hax(1),'xy');  % This means image will be upside-down from how Matlab normally displays images
else
  axis(hax(1),'ij'); 
end
set(hax(1), 'DataAspectRatio', [1 1 1]) ;
set(hax(1), 'XLim', [0.5 size(bkgdImage,2)+0.5]) ;
set(hax(1), 'YLim', [0.5 size(bkgdImage,1)+0.5]) ;

hold(hax(1),'on');

% Plot template match outline
lo_corner_xy = ledXY - templateShapeXY/2 ;
hi_corner_xy = ledXY + templateShapeXY/2 ;
corners = [ lo_corner_xy(1) lo_corner_xy(2) ; ...
            hi_corner_xy(1) lo_corner_xy(2) ; ...
            hi_corner_xy(1) hi_corner_xy(2) ; ...
            lo_corner_xy(1) hi_corner_xy(2) ; ...
            lo_corner_xy(1) lo_corner_xy(2) ] ; 
x_data = corners(:,1)' ;
y_data = corners(:,2)' ;
line('Parent', hax(1), 'XData', x_data, 'YData', y_data, 'Color', 'w', 'LineWidth', 3) ;
line('Parent', hax(1), 'XData', x_data, 'YData', y_data, 'Color', 'k', 'LineWidth', 1) ;  % 0.5 is the default line width

% plot LED (cluster) location
line('Parent', hax(1), ...
     'XData', ledXY(1), ...
     'YData', ledXY(2), ...
     'Marker', 's', ...
     'MarkerSize', 12, ...
     'MarkerEdgeColor', 'k', ...
     'MarkerFaceColor', 'w');
line('Parent', hax(1), ...
     'XData', ledXY(1), ...
     'YData', ledXY(2), ...
     'Marker', 'x', ...
     'MarkerSize', 12, ...
     'MarkerEdgeColor', 'k', ...
     'MarkerFaceColor', 'none');

title(hax(1),'Maxvalue image and LED detection');
colormap(hfig,'jet');

try
  if exist(imsavename,'file'),
    try
      delete(imsavename);
    catch ME,
      warning('Could not delete file %s:\n %s',imsavename,getReport(ME));
    end
  end
  set(hfig,'Units','pixels','Position',figpos);
  save2png(imsavename,hfig);
catch ME,
  warning('Could not save to %s:\n%s',imsavename,getReport(ME));
end
