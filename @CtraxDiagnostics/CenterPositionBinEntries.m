function [heatmap,hfig] = CenterPositionBinEntries(obj,varargin)

% parse inputs
nbins_x = [];
nbins_y = [];
nbins = 100;
xlim = nan(1,2);
ylim = nan(1,2);
doplot = false;
hfig = [];
[nbins_x,nbins_y,nbins,xlim,ylim,doplot,hfig] = ...
  myparse(varargin,'nbins_x',nbins_x,'nbins_y',nbins_y,'nbins',nbins,...
  'xlim',xlim,'ylim',ylim,'doplot',doplot,'hfig',hfig);

% set edges if not input
if isempty(nbins_x),
  nbins_x = nbins;
end
if isnan(xlim(1)),
  xlim(1) = min([obj.trx.x_mm]);
end
if isnan(xlim(2)),
  xlim(2) = max([obj.trx.x_mm]);
end
binwidth = (xlim(2)-xlim(1))/nbins_x;
edges_x = linspace(xlim(1),xlim(2),nbins_x+1);

if isempty(nbins_y),
  nbins_y = nbins;
end
if isnan(ylim(1)),
  ylim(1) = min([obj.trx.y_mm]);
end
if isnan(ylim(2)),
  ylim(2) = max([obj.trx.y_mm]);
end
binheight = (ylim(2)-ylim(1))/nbins_y;
edges_y = linspace(ylim(1),ylim(2),nbins_y+1);

centers_x = (edges_x(1:end-1)+edges_x(2:end))/2;
centers_y = (edges_y(1:end-1)+edges_y(2:end))/2;

nflies = length(obj.trx);

% allocate
heatmap.fracperfly = nan(nbins_y,nbins_x,nflies);
heatmap.counts = zeros(nbins_y,nbins_x);
% number of data points per fly
heatmap.n = nan(1,nflies);

% loop over flies
for fly = 1:nflies,
  
  % find which bin the fly is in
  binx = floor((obj.trx(fly).x_mm - xlim(1))/binwidth);
  biny = floor((obj.trx(fly).y_mm - ylim(1))/binheight);
  
  % find changes in bin
  isnewbin = [true,(binx(1:end-1)~=binx(2:end)) | (biny(1:end-1)~=biny(2:end))];

  % rehistogram -- probably inefficient, but whatever
  countscurr = hist3([obj.trx(fly).y_mm(isnewbin);obj.trx(fly).x_mm(isnewbin)]',{edges_y,edges_x});
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

% store the centers and edges
heatmap.centers_x = centers_x;
heatmap.centers_y = centers_y;
heatmap.edges_x = edges_x;
heatmap.edges_y = edges_y;

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
    set(hfig,'Position',[21,1,1000,800]);
  end
  hax = createsubplots(2,2,.05);
  
  % heatmap of all flies' positions combined
  axes(hax(1)); %#ok<*MAXES>
  imagesc(xlim,ylim,heatmap.fracallflies);
  axis image;
  title('Fraction of bin entries over all flies combined');
  colorbar;

  % heatmap of average fraction per fly
  axes(hax(2));
  imagesc(xlim,ylim,heatmap.meanfracperfly);
  axis image;
  title('Mean fraction of bin entries per fly');
  colorbar;

  % heatmap of all flies' positions combined
  axes(hax(3));
  imagesc(xlim,ylim,log(heatmap.fracallflies));
  axis image;
  title('Log-fraction of bin entries over all flies combined');
  colorbar;

  % heatmap of average fraction per fly
  axes(hax(4));
  imagesc(xlim,ylim,log(heatmap.meanfracperfly));
  axis image;
  title('Log-mean fraction of bin entries per fly');
  colorbar;
  
  linkaxes(hax);

end