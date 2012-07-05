function [centers,edges,meanfrac,stdfrac,stderrfrac,nexpsanalyzed,Zpertype] = ...
  SelectHistData2Exp(fns,bininfo,hist_plot_params,plottype,...
  histperexp,meanhistperexp)

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
nexpsanalyzed = nan(1,ntypes);
Zpertype = nan(1,ntypes);

if nbinscombine == 1,

  for typei = 1:ntypes,
    fn = fns{typei};
    if strcmpi(plottype,'linear'),
      meanfrac{typei} = meanhistperexp.(fn).meanfrac_linear;
      stdfrac{typei} = meanhistperexp.(fn).stdfrac_linear;
    else
      meanfrac{typei} = meanhistperexp.(fn).meanfrac_log;
      stdfrac{typei} = meanhistperexp.(fn).stdfrac_log;
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

  % recompute mean, std from per-exp histograms
  for typei = 1:ntypes,

    fn = fns{typei};
    
    % which exps we will look at
    goodidx = ~isnan(histperexp.(fn).Z) & histperexp.(fn).Z > 0;

    % per-exp hists
    if strcmpi(plottype,'linear'),
      fracperexp0 = histperexp.(fn).meanfrac_linear';
    else
      fracperexp0 = histperexp.(fn).meanfrac_log';
    end
    nexps = size(fracperexp0,2);
    
    % combine bins
    fracperexp = nan(nbins,nexps);
    for i0 = 1:nbins,
      j0 = (i0-1)*nbinscombine+1;
      j1 = min(nbins0,j0+nbinscombine-1);
      fracperexp(i0,:) = sum(fracperexp0(j0:j1,:),1);
    end
    
    % take stats over exps
    meanfrac{typei} = nan(1,nbins);
    stdfrac{typei} = nan(1,nbins);
    for j = 1:nbins,
      [meanfrac{typei}(j)] = ...
        mean(fracperexp(j,goodidx),2);
%       [meanfrac{typei}(j),s] = ...
%         weighted_mean_cov(fracperfly(j,goodidx)',...
%         histperexp.(fn).fracframesanalyzed(goodidx)');
      stdfrac{typei}(j) = std(fracperexp(j,goodidx),1,2);
    end
  end
end
      
% these computations are the same whether we combine bins or not
for typei = 1:ntypes,
  fn = fns{typei};
  goodidx = ~isnan(histperexp.(fn).Z) & histperexp.(fn).Z > 0;
  nexps = nnz(goodidx);
  stderrfrac{typei} = stdfrac{typei} / sqrt(nexps);
  nexpsanalyzed(typei) = meanhistperexp.(fn).n;
  Zpertype(typei) = meanhistperexp.(fn).Z;
end
