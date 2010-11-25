function [heatmap,hfig] = CenterPositionBinEntriesPolar(obj,varargin)

% parse inputs
nbins_r = [];
nbins_theta = [];
nbins = 25;
rlim = [0,nan];
thetalim = [-pi,pi];
doplot = false;
hfig = [];
arena_center = [nan,nan];
xlim = nan(1,2);
ylim = nan(1,2);
[nbins_r,nbins_theta,nbins,rlim,thetalim,xlim,ylim,doplot,hfig,arena_center] = ...
  myparse(varargin,'nbins_r',nbins_r,'nbins_theta',nbins_theta,'nbins',nbins,...
  'rlim',rlim,'thetalim',thetalim,'xlim',xlim,'ylim',ylim,...
  'doplot',doplot,'hfig',hfig,...
  'arena_center',arena_center);

nflies = length(obj.trx);

% set arena center
if isnan(arena_center(1)),
  if isnan(xlim(1)),
    xlim(1) = min([trx.x_mm]);
  end
  if isnan(xlim(2)),
    xlim(2) = max([trx.x_mm]);
  end
  arena_center(1) = mean(xlim);
end
if isnan(arena_center(2)),
  if isnan(ylim(1)),
    ylim(1) = min([trx.y_mm]);
  end
  if isnan(ylim(2)),
    ylim(2) = max([trx.y_mm]);
  end
  arena_center(2) = mean(ylim);
end

% convert to polar coordinates
r = cell(1,nflies);
theta = cell(1,nflies);
for fly = 1:nflies,
  dx = trx(fly).x_mm-arena_center(1);
  dy = trx(fly).y_mm-arena_center(2);
  r{fly} = sqrt( dx.^2 + dy.^2 );
  theta{fly} = atan2(dy,dx);
end

% set edges if not input
if isempty(nbins_r),
  nbins_r = nbins;
end
if isnan(rlim(1)),
  rlim(1) = min([r{:}]);
end
if isnan(rlim(2)),
  rlim(2) = max([r{:}]);
end
binsize_r = (rlim(2)-rlim(1))/nbins_r;
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
heatmap.n = nan(1,nflies);

% loop over flies
for fly = 1:nflies,
  
  % find which bin the fly is in
  binr = floor((r{fly} - rlim(1))/binsize_r);
  bintheta = floor((theta{fly} - thetalim(1))/binsize_theta);
  
  % find changes in bin
  isnewbin = [true,(binr(1:end-1)~=binr(2:end)) | (bintheta(1:end-1)~=bintheta(2:end))];

  % rehistogram -- probably inefficient, but whatever
  countscurr = hist3([r{fly}(isnewbin);theta{fly}(isnewbin)]',{edges_r,edges_theta});
  countscurr(:,end-1) = countscurr(:,end-1)+countscurr(:,end);
  countscurr(end-1,:) = countscurr(end-1,:)+countscurr(end,:);
  countscurr = countscurr(1:end-1,1:end-1);
  
  % normalize per fly
  heatmap.n(fly) = nnz(isnewbin);
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
  if isempty(hfig),
    hfig = figure;
    set(hfig,'position',[21 47 1000 825]);
  else
    figure(hfig);
    clf;
  end
  hax = createsubplots(2,2,.05);
  
  % heatmap of all flies' positions combined
  axes(hax(1)); %#ok<*MAXES>
  polarimagesc(edges_r,edges_theta,arena_center,heatmap.fracallfliespermm2);
  axis image;
  title('Fraction of bin entries over all flies combined per mm^2');
  colorbar;

  % heatmap of average fraction per fly
  axes(hax(2));
  polarimagesc(edges_r,edges_theta,arena_center,heatmap.meanfracperflypermm2);
  axis image;
  title('Mean fraction of bin entries per fly per mm^2');
  colorbar;

  % heatmap of all flies' positions combined
  axes(hax(3));
  polarimagesc(edges_r,edges_theta,arena_center,log(heatmap.fracallfliespermm2));
  axis image;
  title('Log-fraction of bin entries over all flies combined per mm^2');
  colorbar;

  % heatmap of average fraction per fly
  axes(hax(4));
  polarimagesc(edges_r,edges_theta,arena_center,log(heatmap.meanfracperflypermm2));
  axis image;
  title('Log-mean fraction of bin entries per fly per mm^2');
  colorbar;
  
  linkaxes(hax);

end