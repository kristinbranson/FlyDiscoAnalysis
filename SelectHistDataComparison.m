function [centers,edges,meanfrac,stdfrac,stderrfrac,nfliesanalyzed,Zperdata] = ...
  SelectHistDataComparison(fn,bininfo,hist_plot_params,plottype,...
  histperfly,histperexp)

ndata = numel(histperfly);

if strcmpi(plottype,'linear'),
  centers = bininfo.centers_linear;
  edges = bininfo.edges_linear;
else
  centers = bininfo.centers_log;
  edges = bininfo.edges_log;
end

if isfield(hist_plot_params,'nbinscombine'),
  nbinscombine = hist_plot_params.nbinscombine;
else
  nbinscombine = 1;
end

meanfrac = cell(ndata,1);
stdfrac = cell(ndata,1);
stderrfrac = cell(ndata,1);
nfliesanalyzed = nan(ndata,1);
Zperdata = nan(ndata,1);

if nbinscombine == 1,

  for datai = 1:ndata,
    if strcmpi(plottype,'linear'),
      meanfrac{datai} = histperexp{datai}.(fn).meanfrac_linear;
      stdfrac{datai} = histperexp{datai}.(fn).stdfrac_linear;
    else
      meanfrac{datai} = histperexp{datai}.(fn).meanfrac_log;
      stdfrac{datai} = histperexp{datai}.(fn).stdfrac_log;
    end
  end

else
  
  % combine edges, centers
  nbins0 = numel(centers);
  tmp = edges(1:nbinscombine:end);
  if tmp(end) < edges(end),
    tmp(end+1) = edges(end);
  end
  edges = tmp;
  centers = (edges(1:end-1)+edges(2:end))/2;
  nbins = numel(centers);

  % recompute mean, std from per-fly histograms
  for datai = 1:ndata,

    % which flies we will look at
    goodidx = ~isnan(histperfly{datai}.(fn).Z);

    % per-fly hists
    if strcmpi(plottype,'linear'),
      fracperfly0 = histperfly{datai}.(fn).frac_linear;
    else
      fracperfly0 = histperfly{datai}.(fn).frac_log;
    end
    nflies = size(fracperfly0,2);
    
    % combine bins
    fracperfly = nan(nbins,nflies);
    for i0 = 1:nbins,
      j0 = (i0-1)*nbinscombine+1;
      j1 = min(nbins0,j0+nbinscombine-1);
      fracperfly(i0,:) = sum(fracperfly0(j0:j1,:),1);
    end
    
    % take stats over flies
    meanfrac{datai} = nan(1,nbins);
    stdfrac{datai} = nan(1,nbins);
    for j = 1:nbins,
      [meanfrac{datai}(j),s] = ...
        weighted_mean_cov(fracperfly(j,goodidx)',...
        histperfly{datai}.(fn).fracframesanalyzed(goodidx)');
      stdfrac{datai}(j) = sqrt(s);
    end
  end
end
      
% these computations are the same whether we combine bins or not
for datai = 1:ndata,
  goodidx = ~isnan(histperfly{datai}.(fn).Z);
  nflies = sum(histperfly{datai}.(fn).fracframesanalyzed(goodidx));
  stderrfrac{datai} = stdfrac{datai} / sqrt(nflies);
  % backwards compatible
  if isfield(histperexp{datai}.(fn),'nflies') && ~isfield(histperexp{datai}.(fn),'n')
    histperexp{datai}.(fn).n = histperexp{datai}.(fn).nflies;
  end
  nfliesanalyzed(datai) = histperexp{datai}.(fn).n;
  Zperdata(datai) = histperexp{datai}.(fn).Z;
end
