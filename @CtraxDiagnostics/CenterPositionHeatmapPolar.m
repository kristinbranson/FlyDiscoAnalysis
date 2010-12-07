function heatmap = CenterPositionHeatmapPolar(obj,varargin)

% parse inputs
nbins_r = [];
nbins_theta = [];
nbins = 25;
rlim = [0,nan];
thetalim = [-pi,pi];
doplot = false;
hfig = [];
arena_center = obj.arena_center;
arena_radius = obj.arena_radius;
[nbins_r,nbins_theta,nbins,rlim,thetalim,doplot,hfig] = ...
  myparse(varargin,'nbins_r',nbins_r,'nbins_theta',nbins_theta,'nbins',nbins,...
  'rlim',rlim,'thetalim',thetalim,...
  'doplot',doplot,'hfig',hfig);

nflies = length(obj.trx);

% convert to polar coordinates
r = cell(1,nflies);
theta = cell(1,nflies);
for fly = 1:nflies,
  dx = obj.trx(fly).x_mm-arena_center(1);
  dy = obj.trx(fly).y_mm-arena_center(2);
  r{fly} = sqrt( dx.^2 + dy.^2 );
  theta{fly} = atan2(dy,dx);
end

% set edges if not input
if isempty(nbins_r),
  nbins_r = nbins;
end
if isnan(rlim(1)),
  rlim(1) = 0;
end
if isnan(rlim(2)),
  rlim(2) = arena_radius;
end
%binsize_r = (rlim(2)-rlim(1))/nbins_r;
edges_r = linspace(rlim(1),rlim(2),nbins_r+1);

if isempty(nbins_theta),
  nbins_theta = nbins;
end
if isnan(thetalim(1)),
  thetalim(1) = min([theta{:}]);
end
if isnan(thetalim(2)),
  thetalim(2) = max([theta{:}]);
end
binsize_theta = (thetalim(2)-thetalim(1))/nbins_theta;
edges_theta = linspace(thetalim(1),thetalim(2),nbins_theta+1);

centers_r = (edges_r(1:end-1)+edges_r(2:end))/2;
centers_theta = (edges_theta(1:end-1)+edges_theta(2:end))/2;


% allocate
heatmap.fracperfly = nan(nbins_r,nbins_theta,nflies);
heatmap.counts = zeros(nbins_r,nbins_theta);
% number of data points per fly
heatmap.n = [obj.trx.nframes];

% loop over flies
for fly = 1:nflies,
  
  % histogram positions for the current fly
  countscurr = hist3([r{fly};theta{fly}]',{edges_r,edges_theta});
  countscurr(:,end-1) = countscurr(:,end-1)+countscurr(:,end);
  countscurr(end-1,:) = countscurr(end-1,:)+countscurr(end,:);
  countscurr = countscurr(1:end-1,1:end-1);
  
  % normalize per fly
  heatmap.fracperfly(:,:,fly) = countscurr / heatmap.n(fly);
  
  % add to counts
  heatmap.counts = heatmap.counts + countscurr;
  
end

% compute mean, std over flies
heatmap.meanfracperfly = nanmean(heatmap.fracperfly,3);
heatmap.stdfracperfly = nanstd(heatmap.fracperfly,1,3);

% normalize counts over all flies
heatmap.fracallflies = heatmap.counts / sum(heatmap.n);

% area of a bin is (r2^2-r1^2)*binsize_theta/2
binarea_r = (edges_r(2:end).^2-edges_r(1:end-1).^2)*binsize_theta/2;
binarea_r = binarea_r(:);

% normalize by bin area
heatmap.fracperflypermm2 = bsxfun(@rdivide,heatmap.fracperfly,binarea_r);
heatmap.meanfracperflypermm2 = bsxfun(@rdivide,heatmap.meanfracperfly,binarea_r);
heatmap.stdfracperflypermm2 = bsxfun(@rdivide,heatmap.stdfracperfly,binarea_r);
heatmap.fracallfliespermm2 = bsxfun(@rdivide,heatmap.fracallflies,binarea_r);

% store the centers and edges
heatmap.centers_r = centers_r;
heatmap.centers_theta = centers_theta;
heatmap.edges_r = edges_r;
heatmap.edges_theta = edges_theta;

if doplot,
  
  % create the figure and subplots
  doresize = isempty(hfig) || ~ishandle(hfig);
  if isempty(hfig),
    hfig = figure;
  else
    figure(hfig);
    clf;
  end
  if doresize,
    set(hfig,'Position',[4,1,1000,800]);
  end
  heatmap.hfig = hfig;
  hax = createsubplots(2,2,.05);
  heatmap.hax = hax;
  
  % heatmap of all flies' positions combined
  axes(hax(1)); %#ok<*MAXES>
  polarimagesc(edges_r,edges_theta,arena_center,heatmap.fracallfliespermm2);
  axis image;
  title('Fraction of time spent in each bin for all flies combined per mm^2');
  colorbar;

  % heatmap of average fraction per fly
  axes(hax(2));
  polarimagesc(edges_r,edges_theta,arena_center,heatmap.meanfracperflypermm2);
  axis image;
  title('Mean fraction of time spent in each bin per fly per mm^2');
  colorbar;

  % heatmap of all flies' positions combined
  axes(hax(3));
  polarimagesc(edges_r,edges_theta,arena_center,log(heatmap.fracallfliespermm2));
  axis image;
  title('Log-fraction of time spent in each bin for all flies combined per mm^2');
  colorbar;

  % heatmap of average fraction per fly
  axes(hax(4));
  polarimagesc(edges_r,edges_theta,arena_center,log(heatmap.meanfracperflypermm2));
  axis image;
  title('Log-mean fraction of time spent in each bin for per fly per mm^2');
  colorbar;
  
  linkaxes(hax);

end