function saveLedDetectionImage(imsavename, bkgdImage, ledXY, plotYAxisPointsUp)

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
if plotYAxisPointsUp ,
  axis(hax(1),'xy');  % This means image will be upside-down from how Matlab normally displays images
else
  axis(hax(1),'ij'); 
end

hold(hax(1),'on');

% plot bowl marker(s)
plot(hax(1),ledXY(1),ledXY(2),'ks','markersize',12);
plot(hax(1),ledXY(1),ledXY(2),'kx','markersize',12);

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
