function [xfig,yfig] = plot2normfigurecoords(xplot,yplot,hax)

if nargin < 3,
  hax = gca;
end

oldunits = get(hax,'Units');
set(hax,'Units','normalized');
axpos = get(hax,'Position');
if ~strcmpi(oldunits,'normalized'),
  set(hax,'Units',oldunits);
end

xlim = get(hax,'XLim');
ylim = get(hax,'YLim');

dx0plot = xplot - xlim(1);
dy0plot = yplot - ylim(2);

cx = axpos(3)/diff(xlim);
cy = axpos(4)/diff(ylim);

xfig = axpos(1) + dx0plot * cx; 
yfig = axpos(2) + dy0plot * cy; 