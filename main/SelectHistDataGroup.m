function [centers,edges,meanfrac,stdfrac,stderrfrac,nfliesanalyzed,Zpertype] = ...
  SelectHistDataGroup(fns,fields,bins,hist_plot_params,plottype,...
  histperfly,histperexp)

ntypes = numel(fns);

centers = cell(1,ntypes);
edges = cell(1,ntypes);

for i = 1:ntypes,
  centers{i} = bins.(fields{i}).(['centers_',plottype]);
  edges{i} = bins.(fields{i}).(['edges_',plottype]);
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
Zpertype = nan(1,ntypes);

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
  nbins = nan(1,ntypes);
  for typei = 1:ntypes,
    nbins0 = numel(centers{typei});
    tmp = edges{typei}(1:nbinscombine:end);
    if tmp(end) < edges{typei}(end),
      tmp(end+1) = edges{typei}(end); %#ok<AGROW>
    end
    edges{typei} = tmp;
    centers{typei} = (edges{typei}(1:end-1)+edges{typei}(2:end))/2;
    nbins(typei) = numel(centers{typei});

    % recompute mean, std from per-fly histograms
    fn = fns{typei};
    
    % which flies we will look at
    goodidx = ~isnan(histperfly.(fn).Z) & histperfly.(fn).Z > 0;

    % per-fly hists
    if strcmpi(plottype,'linear'),
      fracperfly0 = histperfly.(fn).frac_linear;
    else
      fracperfly0 = histperfly.(fn).frac_log;
    end
    nflies = size(fracperfly0,2);
    
    % combine bins
    fracperfly = nan(nbins(typei),nflies);
    for i0 = 1:nbins(typei),
      j0 = (i0-1)*nbinscombine+1;
      j1 = min(nbins0,j0+nbinscombine-1);
      fracperfly(i0,:) = sum(fracperfly0(j0:j1,:),1);
    end
    
    % take stats over flies
    meanfrac{typei} = nan(1,nbins(typei));
    stdfrac{typei} = nan(1,nbins(typei));
    for j = 1:nbins(typei),
      [meanfrac{typei}(j)] = ...
        weighted_mean(fracperfly(j,goodidx)',...
        histperfly.(fn).Zfly(goodidx)');
      [s] = weighted_scatter(fracperfly(j,goodidx)',...
        meanfrac{typei}(j),histperfly.(fn).fracframesanalyzed(goodidx)');
%       [meanfrac{typei}(j),s] = ...
%         weighted_mean_cov(fracperfly(j,goodidx)',...
%         histperfly.(fn).fracframesanalyzed(goodidx)');
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
  % backwards compatible
  if isfield(histperexp.(fn),'nflies') && ~isfield(histperexp.(fn),'n')
    histperexp.(fn).n = histperexp.(fn).nflies;
  end
  nfliesanalyzed(typei) = histperexp.(fn).n;
  Zpertype(typei) = histperexp.(fn).Z;
end
