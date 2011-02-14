function PlotPerFrameStats(stats_perframefeatures,statsperfly,statsperexp,varargin)

% parse plotting parameters
[hfig,visible,position,axposition,stdalpha] = ...
  myparse(varargin,'hfig',3,...
  'visible','off',...
  'position',[1 1 1000 500],...
  'axposition',[.1,.1,.85,.85],...
  'stdalpha',.2);

% set up figure
if ~ishandle(hfig),
  figure(hfig);
else
  clf(hfig);
end
set(hfig,'Visible',visible,'Position',position);
hax = axes('Position',axposition,'Parent',hfigs);

