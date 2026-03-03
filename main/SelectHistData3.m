function [centers,edges,meanfrac,stdfrac,stderrfrac,nanalyzed,Zpertype] = ...
  SelectHistData3(fns,bininfo,hist_plot_params,plottype,...
  histsingle,histcombined,varargin)

[stattype,meanweighttype,stdweighttype,nframestotal] = myparse(varargin,...
  'stattype','flymeans',...
  'meanweighttype','nframesfly',...
  'stdweighttype','fracframesfly',...
  'nframestotal',[]);

iscelldata = iscell(histsingle);
ndata = numel(histsingle);
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

meanfrac = cell(ndata,ntypes);
stdfrac = cell(ndata,ntypes);
stderrfrac = cell(ndata,ntypes);
nanalyzed = nan(ndata,ntypes);
Zpertype = nan(ndata,ntypes);

if nbinscombine == 1,

  for typei = 1:ntypes,
    fn = fns{typei};
    for datai = 1:ndata,
      if iscelldata,
        datacurr = histcombined{datai}.(fn);
      else
        datacurr = histcombined(datai).(fn);
      end
      if strcmpi(plottype,'linear'),
        meanfrac{datai,typei} = datacurr.meanfrac_linear;
        stdfrac{datai,typei} = datacurr.stdfrac_linear;
      else
        meanfrac{datai,typei} = datacurr.meanfrac_log;
        stdfrac{datai,typei} = datacurr.stdfrac_log;
      end
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

    for datai = 1:ndata,
      if iscelldata,
        datacurr = histsingle{datai}.(fn);
      else
        datacurr = histsingle(datai).(fn);
      end
    
      % which flies we will look at
      goodidx = ~isnan(datacurr.Z) & datacurr.Z > 0;

      % per-fly hists
      if strcmp(stattype,'flymeans'),
        if strcmpi(plottype,'linear'),
          frac0 = datacurr.frac_linear;
        else
          frac0 = datacurr.frac_log;
        end
      else
        if strcmpi(plottype,'linear'),
          frac0 = datacurr.meanfrac_linear';
        else
          frac0 = datacurr.meanfrac_log';
        end
      end
      nsingle = size(frac0,2);
    
      % combine bins
      frac = nan(nbins,nsingle);
      for i0 = 1:nbins,
        j0 = (i0-1)*nbinscombine+1;
        j1 = min(nbins0,j0+nbinscombine-1);
        frac(i0,:) = sum(frac0(j0:j1,:),1);
      end
    
      % take stats over flies
      meanfrac{datai,typei} = nan(1,nbins);
      stdfrac{datai,typei} = nan(1,nbins);
      if strcmpi(stattype,'flymeans'),
        switch lower(meanweighttype),
          case 'nframesfly'
            wmean = datacurr.Zfly(goodidx);
          case 'fracframesfly',
            wmean = datacurr.fracframesanalyzed(goodidx);
          case 'nframesanalyzed',
            wmean = datacurr.Z(goodidx);
          case 'uniform',
            wmean = ones(1,nnz(goodidx));
        end
      else
        wmean = ones(1,nnz(goodidx));
      end
      if strcmpi(stattype,'flymeans'),
        switch lower(stdweighttype),
          case 'nframesfly'
            wstd = datacurr.Zfly(goodidx);
          case 'fracframesfly',
            wstd = datacurr.fracframesanalyzed(goodidx);
          case 'nframesanalyzed',
            wstd = datacurr.Z(goodidx);
          case 'uniform',
            wstd = ones(1,nnz(goodidx));
        end
      else
        wstd = ones(1,nnz(goodidx));
      end
      for j = 1:nbins,
        [meanfrac{datai,typei}(j)] = ...
          weighted_mean(frac(j,goodidx)',wmean');
        [s] = weighted_scatter(frac(j,goodidx)',meanfrac{datai,typei}(j),wstd');
        stdfrac{datai,typei}(j) = sqrt(s);
      end
    end
  end
end
      
% these computations are the same whether we combine bins or not
for typei = 1:ntypes,
  fn = fns{typei};
  for datai = 1:ndata,
    if iscelldata,
      datacurr = histsingle{datai}.(fn);
    else
      datacurr = histsingle(datai).(fn);
    end
    goodidx = ~isnan(datacurr.Z);
    
    if strcmpi(stattype,'flymeans'),
      switch lower(stdweighttype),
        case 'nframesfly'
          if isempty(nframestotal),
            error('nframestotal required for weighttype %s',weighttype);
          end
          n = sum(datacurr.Zfly(goodidx)) / nframestotal;
        case 'fracframesfly',
          n = sum(datacurr.fracframesanalyzed(goodidx));
        case 'nframesanalyzed',
          n = nan;
        case 'uniform',
          n = nnz(goodidx);
      end
    else
      n = nnz(goodidx);
    end
    stderrfrac{datai,typei} = stdfrac{datai,typei} / sqrt(n);
    if iscelldata,
      datacurr = histcombined{datai}.(fn);
    else
      datacurr = histcombined(datai).(fn);
    end
    nanalyzed(datai,typei) = datacurr.n;
    Zpertype(datai,typei) = datacurr.Z;
  end
end
