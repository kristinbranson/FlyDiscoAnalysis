function biasmap = BowlAngleBias(obj,varargin)

% parse inputs
doplot = false;
hfig = [];
arena_center = [0,0];
nbins = 360;
anglerange = pi;
[doplot,hfig,arena_center,nbins,anglerange] = ...
  myparse(varargin,'doplot',doplot,'hfig',hfig,...
  'arena_center',arena_center,'nbins',nbins,...
  'anglerange',anglerange);

nflies = length(obj.trx);

% polar angle
theta = cell(1,nflies);
for fly = 1:nflies,
  dx = obj.trx(fly).x_mm-arena_center(1);
  dy = obj.trx(fly).y_mm-arena_center(2);
  %r{fly} = sqrt( dx.^2 + dy.^2 );
  theta{fly} = atan2(dy,dx);
end

% bins
bin_edges = linspace(-pi,pi,nbins+1);
bin_width = bin_edges(2)-bin_edges(1);
bin_centers = (bin_edges(1:end-1)+bin_edges(2:end))/2;

% how many bins fall within range
nbinsin = floor(anglerange / bin_width);
nbinsout = nbins - nbinsin;
fil = ones(1,nbinsin);

% normalize so that the areas in the two sides are the same
Z = nbinsin / nbinsout;

% center of the range
% nwithin(nbinsin) = sum(counts(1:nbinsin))
% so thetares(nbinsin) = (theta(1)+theta(nbinsin))/2
% = (theta(nbinsin) - bin_width*(nbinsin-1) + theta(nbinsin)) /2
% = theta(nbinsin) - bin_width*(nbinsin-1)/2
biasmap.theta = modrange(bin_centers - bin_width*(nbinsin-1)/2,-pi,pi);
% total number within
nwithinallflies = zeros(1,nbins);
nwithinperfly = zeros(1,nbins,nflies);
% number within per fly
biasmap.biasperfly = zeros([1,nbins,nflies]);
% number of frames per fly
biasmap.n = [obj.trx.nframes];

for fly = 1:nflies,

  n = length(theta{fly});
  
  % histogram into a large number of bins
  counts = hist(theta{fly},bin_centers);
  
  % circular filter to get sum within each sequence of nbinsin bins
  % nwithin(nbinsin) = sum(counts(1:nbinsin))
  nwithin = cconv(fil,counts,nbins);
  nwithout = n - nwithin;
  
  % normalize for this fly
  biasmap.biasperfly(:,:,fly) = (nwithin - nwithout*Z)/n;
  
  % add into counts for all flies
  nwithinallflies = nwithinallflies + nwithin;
  nwithinperfly(:,:,fly) = nwithin;

end

nallflies = sum(biasmap.n);
nwithoutallflies = nallflies - nwithinallflies;
biasmap.biasallflies = (nwithinallflies - nwithoutallflies*Z)/nallflies;
biasmap.biasallflies_jackknife_std = zeros(1,nbins);
for i = 1:nflies,

  idx = true(1,nflies); idx(i) = false;
  nwithinallflies_holdout = sum(nwithinperfly(:,:,idx),3);
  n = sum(biasmap.n(idx));
  nwithoutallflies_holdout = n - nwithinallflies_holdout;
  biasallflies_holdout = (nwithinallflies_holdout - nwithoutallflies_holdout*Z)/nallflies;
  biasmap.biasallflies_jackknife_std = biasmap.biasallflies_jackknife_std + ...
    (biasallflies_holdout - biasmap.biasallflies).^2;
  
end
biasmap.biasallflies_jackknife_std = sqrt(biasmap.biasallflies_jackknife_std/nflies);
biasmap.meanbiasperfly = nanmean(biasmap.biasperfly,3);
biasmap.stdbiasperfly = nanstd(biasmap.biasperfly,1,3);

if doplot,
  doresize = isempty(hfig) || ~ishandle(hfig);
  if isempty(hfig),
    hfig = figure;
  else
    figure(hfig);
    clf;
  end
  if doresize,
    set(hfig,'Position',[101,1,1000,500]);
  end
  
  hax = createsubplots(1,2,.01);
  axes(hax(1)); %#ok<*MAXES>
  polarimagesc([0,min(obj.width_mm,obj.height_mm)],bin_edges,arena_center,biasmap.meanbiasperfly);  
  axis image;
  title(sprintf('Mean bias per fly, anglerange = %.1f deg',anglerange*180/pi));
  colorbar;

  axes(hax(2));
  polarimagesc([0,min(obj.width_mm,obj.height_mm)],bin_edges,arena_center,biasmap.biasallflies);  
  axis image;
  title(sprintf('Bias over all flies, anglerange = %.1f deg',anglerange*180/pi));
  colorbar;

end