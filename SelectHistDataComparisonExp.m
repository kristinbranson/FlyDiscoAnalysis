function [centers,edges,meanfrac,stdfrac,stderrfrac,nexpsanalyzed,Zperdata] = ...
  SelectHistDataComparisonExp(fn,bininfo,hist_plot_params,plottype,...
  histperexp,meanhistperexp)

ndata = numel(histperexp);

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
nexpsanalyzed = nan(ndata,1);
Zperdata = nan(ndata,1);

if nbinscombine == 1,

  for datai = 1:ndata,
    if strcmpi(plottype,'linear'),
      meanfrac{datai} = meanhistperexp{datai}.(fn).meanfrac_linear;
      stdfrac{datai} = meanhistperexp{datai}.(fn).stdfrac_linear;
    else
      meanfrac{datai} = meanhistperexp{datai}.(fn).meanfrac_log;
      stdfrac{datai} = meanhistperexp{datai}.(fn).stdfrac_log;
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
  for datai = 1:ndata,

    % which exps we will look at
    goodidx = ~isnan(histperexp{datai}.(fn).Z) & histperexp{datai}.(fn).Z > 0;

    % per-exp hists
    if strcmpi(plottype,'linear'),
      fracperexp0 = histperexp{datai}.(fn).meanfrac_linear';
    else
      fracperexp0 = histperexp{datai}.(fn).meanfrac_log';
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
    meanfrac{datai} = nan(1,nbins);
    stdfrac{datai} = nan(1,nbins);
    for j = 1:nbins,
      [meanfrac{datai}(j)] = ...
        mean(fracperexp(j,goodidx),2);
      stdfrac{datai}(j) = std(fracperexp(j,goodidx),1,2);
    end
  end
end
      
% these computations are the same whether we combine bins or not
for datai = 1:ndata,
  goodidx = ~isnan(histperexp{datai}.(fn).Z) & histperexp{datai}.(fn).Z > 0;
  nexps = nnz(goodidx);
  stderrfrac{datai} = stdfrac{datai} / sqrt(nexps);
  nexpsanalyzed(datai) = meanhistperexp{datai}.(fn).n;
  Zperdata(datai) = meanhistperexp{datai}.(fn).Z;
end
