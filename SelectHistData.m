function [centers,edges,meanfrac,stdfrac,stderrfrac,nfliesanalyzed] = ...
  SelectHistData(fns,bininfo,hist_plot_params,plottype,...
  histperfly,histperexp)

ntypes = numel(fns);

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

meanfrac = cell(1,ntypes);
stdfrac = cell(1,ntypes);
stderrfrac = cell(1,ntypes);
nfliesanalyzed = nan(1,ntypes);

if nbinscombine == 1,

  for typei = 1:ntypes,
    fn = fns{typei};
    if strcmpi(plottype,'linear'),
      meanfrac{typei} = histperexp.(fn).meanfrac_linear;
      stdfrac{typei} = histperexp.(fn).stdfrac_linear;
    else
      meanfrac{typei} = histperexp.(fn).meanfrac_log;
      stdfrac{typei} = histperexp.(fn).stdfrac_log;
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
  for typei = 1:ntypes,

    fn = fns{typei};
    
    % which flies we will look at
    goodidx = ~isnan(histperfly.(fn).Z);

    % per-fly hists
    if strcmpi(plottype,'linear'),
      fracperfly0 = histperfly.(fn).frac_linear;
    else
      fracperfly0 = histperfly.(fn).frac_log;
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
    meanfrac{typei} = nan(1,nbins);
    stdfrac{typei} = nan(1,nbins);
    for j = 1:nbins,
      [meanfrac{typei}(j),s] = ...
        weighted_mean_cov(fracperfly(j,goodidx)',...
        histperfly.(fn).fracframesanalyzed(goodidx)');
      stdfrac{typei}(j) = sqrt(s);
    end
  end
end
      
% these computations are the same whether we combine bins or not
for typei = 1:ntypes,
  fn = fns{typei};
  goodidx = ~isnan(histperfly.(fn).Z);
  nflies = sum(histperfly.(fn).fracframesanalyzed(goodidx));
  stderrfrac{typei} = stdfrac{typei} / sqrt(nflies);
  nfliesanalyzed(typei) = histperexp.(fn).nflies;
end
