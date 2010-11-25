function heatmap = CenterPositionHeatmap(obj,varargin)

% parse inputs
nbins = 100;
[edges_x,edges_y,nbins_x,nbins_y,nbins,xlim,ylim] = ...
  myparse(varargin,'edges_x',[],'edges_y',[],'nbins_x',[],'nbins_y',[],'nbins',nbins,...
  'xlim',nan(1,2),'ylim',nan(1,2));

% set edges if not input
if isempty(edges_x),
  if isempty(nbins_x),
    nbins_x = nbins;
  end
  if isnan(xlim(1)),
    xlim(1) = min([trx.x_mm]);
  end
  if isnan(xlim(2)),
    xlim(2) = max([trx.x_mm]);
  end
  edges_x = linspace(xlim(1),xlim(2),nbins_x+1);
else
  nbins_x = length(edges_x) - 1;  
end
if isempty(edges_y),
  if isempty(nbins_y),
    nbins_y = nbins;
  end
  if isnan(ylim(1)),
    ylim(1) = min([trx.y_mm]);
  end
  if isnan(ylim(2)),
    ylim(2) = max([trx.y_mm]);
  end
  edges_y = linspace(ylim(1),ylim(2),nbins_y+1);
else
  nbins_y = length(edges_y) - 1;  
end
centers_x = (edges_x(1:end-1)+edges_x(2:end))/2;
centers_y = (edges_y(1:end-1)+edges_y(2:end))/2;

nflies = length(obj.trx);

% allocate
heatmap.fracperfly = nan(nbins_y,nbins_x,nflies);
heatmap.counts = zeros(nbins_y,nbins_x);
% number of data points per fly
heatmap.n = [obj.trx.nframes];

% loop over flies
for fly = 1:nflies,
  
  % histogram positions for the current fly
  countscurr = hist3([obj.trx(fly).y_mm;obj.trx(fly).x_mm]',{edges_y,edges_x});
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

% store the centers and edges
heatmap.centers_x = centers_x;
heatmap.centers_y = centers_y;
heatmap.edges_x = edges_x;
heatmap.edges_y = edges_y;