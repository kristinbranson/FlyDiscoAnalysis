function res = HistogramMeasurement(obj,fn,varargin)

res = struct;

%% parse inputs

if obj.nexpdirs == 0,
  error('No experiments loaded');
end
if ~isfield(obj.trx,fn),
  error('Field %s not computed for trx',fn);
end

[edges,nbins,lim,lim_prctile,...
  hax,hfig,doplot,...
  conditions,outputfun,outputfun_string,...
  expdirs,flies,averaging,jackknife,...
  histstyleparams,perflyhiststyleparams,perexphiststyleparams,errstyleparams,...
  legendstyleparams,axesstyleparams,...
  ploterrorbars,plotperfly,plotperexp,axisplot,docomputestd,docomputestderr,...
  binmode,figpos,ylogscale] = ...
  myparse(varargin,...
  'edges',[],...
  'nbins',obj.histogrammeasurements_nbins,...
  'lim',nan(1,2),...
  'lim_prctile',[1,99],...
  'hax',[],...
  'hfig',[],...
  'doplot',obj.histogrammeasurements_doplot,...
  'conditions',[],...
  'outputfun',[],...
  'outputfun_string','',...
  'expdirs',obj.expdir_bases,...
  'flies',1:obj.nflies,...
  'averaging',obj.histogrammeasurements_averaging,...
  'jackknife',obj.histogrammeasurements_jackknife,...
  'histstyleparams',obj.histogrammeasurements_histstyleparams,...
  'perflyhiststyleparams',obj.histogrammeasurements_perflyhiststyleparams,...
  'perexphiststyleparams',obj.histogrammeasurements_perexphiststyleparams,...
  'errstyleparams',obj.histogrammeasurements_errstyleparams,...
  'legendstyleparams',obj.histogrammeasurements_legendstyleparams,...
  'axesstyleparams',obj.histogrammeasurements_axesstyleparams,...
  'ploterrorbars',obj.histogrammeasurements_ploterrorbars,...
  'plotperfly',obj.histogrammeasurements_plotperfly,...
  'plotperexp',obj.histogrammeasurements_plotperexp,...
  'axisplot',nan(1,4),...
  'docomputestd',obj.histogrammeasurements_docomputestd,...
  'docomputestderr',obj.histogrammeasurements_docomputestderr,...
  'binmode','linear',...
  'figpos',obj.histogrammeasurements_figpos,...
  'ylogscale',true);

% take the intersection of specified flies and expdirs

% which expdirs will we look at?
[didfind,ns] = ismember(expdirs,obj.expdir_bases);
if ~all(didfind),
  warning(['The following expdirs are not loaded: ',sprintf('%s ',expdirs{~didfind})]);
end
ns = ns(didfind);

% which flies does this correspond to?
allflies_perexp = [obj.movie2flies{ns}];
flies = intersect(flies,allflies_perexp);
ns = unique(obj.fly2movie(flies));
nflies = length(flies);
% save for returning
res.ns = ns;
res.nflies = nflies;

% check that there is at least one experiment
nexpdirs = length(ns);
if nexpdirs == 0,
  error('No experiments selected.');
end

% set jackknife method if not set
if isempty(jackknife),
  if nexpdirs > 1,
    jackknife = 'perexp';
  elseif nflies == 1,
    jackknife = 'none';
  else
    jackknife = 'perfly';
  end
end

% to split over experiments, we need more than one experiment
if nexpdirs == 1 && strcmpi(jackknife,'perexp'),
  warning('Only one experiment selected, but jackknife = ''perexp''. Splitting over flies instead.');
  jackknife = 'perfly';
end
% to split over flies, we need more than one fly
if nflies == 1 && strcmpi(jackknife,'perfly'),
  warning('Only one fly selected, but jackknife = ''perfly''. Splitting over frames not implemented. Not jackknifing.');
  jackknife = 'none';
end
res.jackknife = jackknife;
  
% set edges if not input
if isempty(edges),
  if any(isnan(lim)),
    minx = inf;
    maxx = -inf;
    % for now, don't worry about memory, store all the data
    useprctile = any(~isnan(lim_prctile));
    if useprctile,
      xall = [];
    end
    for fly = flies,
      if ~isempty(conditions),
        x = obj.trx(fly).(fn)(conditions(obj.trx(fly)));
      else
        x = obj.trx(fly).(fn);
      end
      if ~isempty(outputfun),
        x = outputfun(x);
      end
      badidx = isinf(x) | isnan(x);
      minx = min(minx,min(x(~badidx)));
      maxx = max(maxx,max(x(~badidx)));
      if useprctile,
        xall = [xall,x]; %#ok<AGROW>
      end
    end
    if isnan(lim(1)),
      if isnan(lim_prctile(1)),
        lim(1) = minx;
      else
        lim(1) = prctile(xall,lim_prctile(1));
      end
    end
    if isnan(lim(2)),
      if isnan(lim_prctile(2)),
        lim(2) = maxx;
      else
        lim(2) = prctile(xall,lim_prctile(2));
      end
    end
  end
  switch lower(binmode),
    case 'linear',
      edges = linspace(lim(1),lim(2),nbins+1);
    case 'log',
      % make sure >= 1
      tmplim = lim - lim(1) + 1;
      edges = exp(linspace(log(tmplim(1)),log(tmplim(2)),nbins+1));
      edges = edges + lim(1) - 1;
    case 'logabs',
      % from lim(1) to 0, we do log spacing on neg value
      % from 0 to lim(2), we do log spacing on pos value
      
      if lim(1) > 0,
        % corner case: both are positive
        tmplim = lim - lim(1) + 1;
        edges = exp(linspace(log(tmplim(1)),log(tmplim(2)),nbins+1));
        edges = edges + lim(1) - 1;

      elseif lim(2) < 0,
        % corner case: both are negative
        tmplim = -lim;
        tmplim = tmplim - tmplim(1) + 1;
        edges = exp(linspace(log(tmplim(1)),log(tmplim(2)),nbins+1));
        edges = fliplr(-(edges - lim(2) - 1));

      else
        % how much of data is below 0
        fracneg = -lim(1) / (lim(2)-lim(1));
        fracpos = 1 - fracneg;
        % how many bins will we have on one side of 0
        if fracneg > .5,
          nbinsneg = floor(fracneg*nbins);
          nbinspos = nbins - nbinsneg;
        else
          nbinspos = floor(fracpos*nbins);
          nbinsneg = nbins - nbinspos;
        end
        % positive edges
        tmplim = [0,lim(2)] + 1;
        edgespos = exp(linspace(log(tmplim(1)),log(tmplim(2)),nbinspos+1));
        edgespos = edgespos - 1;
        % negative edges
        tmplim = [0,-lim(1)];
        tmplim = tmplim + 1;
        edgesneg = exp(linspace(log(tmplim(1)),log(tmplim(2)),nbinsneg+1));
        edgesneg = fliplr(-(edgesneg - 1));
        edges = [edgesneg(1:end-1),edgespos];
      end        
        
    otherwise,
      error('Unknown binmode %s',binmode);
  end
else
  nbins = length(edges) - 1;  
end

centers = (edges(1:end-1)+edges(2:end))/2;

res.edges = edges;
res.centers = centers;
res.nbins = nbins;

%% histogram for each fly
countsperfly = nan(length(flies),nbins);
for i = 1:nflies,
  fly = flies(i);
  if ~isempty(conditions),
    x = obj.trx(fly).(fn)(conditions(obj.trx(fly)));
  else
    x = obj.trx(fly).(fn);
  end
  if ~isempty(outputfun),
    x = outputfun(x);
  end
  countscurr = histc(x,edges);
  countsperfly(i,:) = [countscurr(1:end-2),countscurr(end-1)+countscurr(end)];
end

%% normalize by bin size
% donormbybinwidth = ~ismember(lower(binmode),{'linear'});
% if donormbybinwidth,
%   binwidths = diff(edges);
%   mean_binwidth = mean(binwidths);
%   countsperfly = bsxfun(@rdivide,countsperfly,binwidths/mean_binwidth);
% end

%% set up to index into data structures only computed for selected flies, selected exps

% flyidx i corresponds to fly flies(i)
fly2idx = sparse(ones(1,nflies),flies,1:nflies,1,obj.nflies);
% expidx i corresponds to experiment ns(i)
n2idx = sparse(ones(1,nexpdirs),ns,1:nexpdirs,1,obj.nexpdirs);

expdirs = obj.expdir_bases(ns);
movie2flies = cell(1,nexpdirs);
nfliespermovie = zeros(1,nexpdirs);
for i = 1:length(ns),
  n = ns(i);
  movie2flies{i} = full(fly2idx(obj.movie2flies{n}));
  movie2flies{i}(movie2flies{i} == 0) = [];
  nfliespermovie(i) = length(movie2flies{i});
end
fly2movie = full(n2idx(obj.fly2movie(flies)));

% replace flies, ns with indexes
%flies0 = flies;
%ns0 = ns;
%flies = 1:nflies;
%ns = 1:nexpdirs;

%% per-fly histograms

% these won't all be used for all averaging conditions, but it should be
% quick and will speed up jackknifing, so compute them
Zperfly = sum(countsperfly,2);
fracperfly = bsxfun(@rdivide,countsperfly,Zperfly);

%% compute frac using all the data

% compute histogram using all the data
frac = compute_frac(1:nflies,1:nexpdirs);

res.frac = frac;

%% jackknife to get an estimate of standard deviation

if docomputestd && ~strcmpi(jackknife,'none'),

  jackknife_frac_mean = zeros(1,nbins);
  jackknife_frac_std = zeros(1,nbins);
  nsets = zeros(1,nbins);
  switch jackknife,
    case 'perfly',
      for fly = 1:nflies,
        % use only the current fly
        frac_curr = compute_frac(fly,1:nexpdirs);
        % make sure not nan
        isdata = ~isnan(frac_curr);
        nsets(isdata) = nsets(isdata) + 1;
        % update mean, std estimates
        jackknife_frac_mean(isdata) = jackknife_frac_mean(isdata) + frac_curr(isdata);
        jackknife_frac_std(isdata) = jackknife_frac_std(isdata) + frac_curr(isdata).^2;
      end
    case 'perexp',
      for n = 1:nexpdirs,
        % use only the current exp
        frac_curr = compute_frac(1:nflies,n);
        % make sure not nan
        isdata = ~isnan(frac_curr);
        nsets(isdata) = nsets(isdata) + 1;
        % update mean, stderr estimates
        jackknife_frac_mean(isdata) = jackknife_frac_mean(isdata) + frac_curr(isdata);
        jackknife_frac_std(isdata) = jackknife_frac_std(isdata) + frac_curr(isdata).^2;
      end
    otherwise
      error('Unknown jackknife bagging unit %s',jackknife);
  end
  % normalize
  jackknife_frac_mean = jackknife_frac_mean ./ nsets;
  jackknife_frac_std = sqrt(jackknife_frac_std ./ nsets - jackknife_frac_mean.^2);
end

res.jackknife_frac_std = jackknife_frac_std;

%% jackknife to get an estimate of standard error
if docomputestderr && ~strcmpi(jackknife,'none'),

  jackknife_frac_mean = zeros(1,nbins);
  jackknife_frac_stderr = zeros(1,nbins);
  nsets = zeros(1,nbins);
  switch jackknife,
    case 'perfly',
      for fly = 1:nflies,
        % use all but the current fly
        fliescurr = setdiff(1:nflies,fly);
        frac_curr = compute_frac(fliescurr,1:nexpdirs);
        % make sure not nan
        isdata = ~isnan(frac_curr);
        nsets(isdata) = nsets(isdata) + 1;
        % update mean, stderr estimates
        jackknife_frac_mean(isdata) = jackknife_frac_mean(isdata) + frac_curr(isdata);
        jackknife_frac_stderr(isdata) = jackknife_frac_stderr(isdata) + frac_curr(isdata).^2;
      end
    case 'perexp',
      for n = 1:nexpdirs,
        % use only the current exp
        expscurr = setdiff(1:nexpdirs,n);
        frac_curr = compute_frac(1:nflies,expscurr);
        % make sure not nan
        isdata = ~isnan(frac_curr);
        nsets(isdata) = nsets(isdata) + 1;
        % update mean, stderr estimates
        jackknife_frac_mean(isdata) = jackknife_frac_mean(isdata) + frac_curr(isdata);
        jackknife_frac_stderr(isdata) = jackknife_frac_stderr(isdata) + frac_curr(isdata).^2;
      end
    otherwise
      error('Unknown jackknife bagging unit %s',jackknife);
  end
  % normalize
  jackknife_frac_mean = jackknife_frac_mean ./ nsets;
  jackknife_frac_stderr = jackknife_frac_stderr ./ nsets - jackknife_frac_mean.^2;
  jackknife_frac_stderr = sqrt(jackknife_frac_stderr .* (nsets-1)./nsets);
end

res.jackknife_frac_stderr = jackknife_frac_stderr;
    
%% plot

res.plotparams = struct;

if doplot,
  
  plotinfo = struct;
  plotinfo.expdirs = expdirs;
  plotinfo.flies = flies;
  plotinfo.averaging = averaging;
  plotinfo.jackknife = jackknife;
  
  %% get a handle for axes to plot in
  if isempty(hax),
    if isempty(hfig),
      hfig = figure;
      if ~isempty(figpos),
        set(hfig,'Position',figpos);
      end
    elseif ishandle(hfig),
      clf(hfig);
    else
      figure(hfig);
      if ~isempty(figpos),
        set(hfig,'Position',figpos);
      end
    end
    hax = get(hfig,'CurrentAxes');
    if isempty(hax),
      hax = axes('parent',hfig);
    end
  else
    hfig = get(hax,'Parent');
  end
  cla(hax);
  hold(hax,'on');
  set(hax,'UserData',plotinfo);
  if ylogscale,
    set(hax,'YScale','log');
  else
    set(hax,'YScale','linear');
  end
%   if ismember(lower(binmode),{'linear'}),
%     set(hax,'XScale','linear');
%   else
%     set(hax,'XScale','log');
%   end
  res.plotparams.hax = hax;
  res.plotparams.hfig = hfig;
  
  %% legend stuff
  
  legend_handles = [];
  legend_strings = {};
  
  % which flies, exps will be in the legend
  if nexpdirs > 1,
    legendexps = [1,nexpdirs];
  else
    legendexps = 1;
  end
  legendflies = [];
  for i = 1:length(legendexps),
    legendflies = [legendflies,movie2flies{legendexps(i)}(1)]; %#ok<AGROW>
    if nfliespermovie(legendexps(i)) > 1,
      legendflies = [legendflies,movie2flies{legendexps(i)}(end)]; %#ok<AGROW>
    end
  end
  
  %% plot standard error/standard deviation
  
  if ~strcmpi(ploterrorbars,'none') && ~strcmpi(jackknife,'none')
    switch ploterrorbars,
      case 'std',
        err = jackknife_frac_std;
      case 'stderr',
        err = jackknife_frac_stderr;
      otherwise
        error('Unknown ploterrorbars type %s',ploterrorbars);
    end
    
    errinfo = struct;
    errinfo.jackknife = jackknife;
    if strcmpi(jackknife,'perfly'),
      errinfo.n = nflies;
    elseif strcmpi(jackknife,'perexp'),
      errinfo.n = nexpdirs;
    end
    
    herr = patch([centers,fliplr(centers)],[frac-err,fliplr(frac+err)],.8*ones(1,3),...
      'parent',hax,'edgecolor','none','tag',sprintf('err_%s',ploterrorbars),'userdata',errinfo);
    if ~isempty(errstyleparams),
      set(herr,errstyleparams{:});
    end
    
    legend_handles(end+1) = herr; 
    legend_strings{end+1} = sprintf('%s (bagging %s)',ploterrorbars,jackknife);
    
  end

  %% get colors per fly, per experiment
  
  % skip ncolorsskip between experiments
  ncolorsskip = 5;
  % generate colors + extras from colormap
  ncolors = nflies + (nexpdirs-1)*ncolorsskip;
  tmpcolors = jet(ncolors);
  colorsperfly = repmat(tmpcolors(1,:),[nflies,1]);
  i = 1;
  for fly = 2:nflies,
    i = i + 1;
    if fly2movie(fly) ~= fly2movie(fly-1),
      i = i + ncolorsskip;
    end
    colorsperfly(fly,:) = tmpcolors(i,:);
  end
  % color for experiment is the mean of the colors for the corresponding
  % flies
  colorsperexp = zeros(nexpdirs,3);
  for n = 1:nexpdirs,
    colorsperexp(n,:) = mean(colorsperfly(movie2flies{n},:),1);
  end
  % make flies lighter, experiments darker
  colorsperfly = colorsperfly*.8 + .2;
  colorsperexp = colorsperexp*.8;
  res.plotparams.colorsperfly = colorsperfly;
  res.plotparams.colorsperexp = colorsperexp;
  
  % for setting axis limits
  minv = inf;
  maxv = -inf;
  
  %% plot the per-fly histograms
  if plotperfly,
    hperflyhist = zeros(1,nflies);
    for fly = 1:nflies,
      flyglobal = flies(fly);
      nglobal = obj.fly2movie(flyglobal);
      expdir = obj.expdir_bases{nglobal};
      flyofmovie = find(flyglobal == obj.movie2flies{nglobal},1);
      flyinfo = struct('fly',flyglobal,'n',nglobal,'expdir',expdir,'ndata',Zperfly(fly),'flyofmovie',flyofmovie);
      hperflyhist(fly) = plot(hax,centers,fracperfly(fly,:),':','color',colorsperfly(fly,:),'linewidth',.5,'userdata',flyinfo,'tag',sprintf('perfly%d',fly));
    end
    
    % update limits
    minv = min(minv,min(fracperfly(:)));
    maxv = max(maxv,max(fracperfly(:)));

    % set line style
    if ~isempty(perflyhiststyleparams),
      set(hperflyhist,perflyhiststyleparams{:});
    end
    res.plotparams.hperflyhist = hperflyhist;
    
    for fly = legendflies,
      flyglobal = flies(fly);
      nglobal = obj.fly2movie(flyglobal);
      expdir = obj.expdir_bases{nglobal};
      flyofmovie = find(flyglobal == obj.movie2flies{nglobal},1);
      legend_handles(end+1) = hperflyhist(fly); %#ok<AGROW>
      legend_strings{end+1} = sprintf('fly %d of %s',flyofmovie,expdir); %#ok<AGROW>
    end
    
  end
  
  %% plot per-experiment histogram
  
  if plotperexp,
    hperexphist = zeros(1,nexpdirs);
    for n = 1:nexpdirs,
      expinfo = struct('n',ns(n),'flies',flies(movie2flies{n}),'expdir',expdirs{n},'ndataperfly',Zperfly(movie2flies{n}));
      fracperexpcurr = compute_frac(1:nflies,n);
      % update limits
      minv = min(minv,min(fracperexpcurr));
      maxv = max(maxv,max(fracperexpcurr));
      hperexphist(n) = plot(hax,centers,fracperexpcurr,'-','color',colorsperexp(n,:),'linewidth',1,'tag',sprintf('exp%d',n),'userdata',expinfo);
    end
    if ~isempty(perexphiststyleparams),
      set(hperexphist,perexphiststyleparams{:});
    end
    res.plotparams.hperexphist = hperexphist;
    
    for n = legendexps,
      legend_handles(end+1) = hperexphist(n); %#ok<AGROW>
      legend_strings{end+1} = sprintf('exp %s',expdirs{n}); %#ok<AGROW>
    end
    
  end
  
  %% plot the histogram
  hhist = plot(hax,centers,frac,'.-','color',[0,0,0],'linewidth',2,'tag','mainhist','userdata',plotinfo);
  % update limits
  minv = min(minv,min(frac));
  maxv = max(maxv,max(frac));

  if ~isempty(histstyleparams),
    set(hhist,histstyleparams{:});
  end
  res.plotparams.hhist = hhist;
  
  legend_handles(end+1) = hhist;
  legend_strings{end+1} = 'All flies, all experiments';

  %% set axis limits
  if isnan(axisplot(1)),
    axisplot(1) = edges(1);
  end
  if isnan(axisplot(2)),
    axisplot(2) = edges(end);
  end

  dv = maxv - minv;
  if isnan(axisplot(3)),
    axisplot(3) = max(0,minv - dv*.01);
  end
  if isnan(axisplot(4)),
    axisplot(4) = maxv + dv*.01;
  end
  % make sure the y lims are positive
  if ylogscale,
    axisplot(3:4) = max(axisplot(3:4),min(frac(frac>0))/2);
  end
  axis(hax,axisplot);
  
  res.plotparams.axisplot = axisplot;

  %% labels
  
  if ~isempty(outputfun),
    if ~isempty(outputfun_string),
      s = sprintf('%s of %s',outputfun_string,fn);
    else
      s = sprintf('%s of %s',func2str(outputfun),fn);
    end
  else
    s = fn;
  end
  hxlabel = xlabel(hax,s,'interpreter','none');
  set(hfig,'name',s);
  switch averaging,
    case 'allexps_allflies',
      s = 'Fraction of frames';
    case 'allexps_perfly',
      s = 'Mean fraction of frames per fly';
    case 'perexp_allflies',
      s = 'Mean fraction of frames per experiment';
    case 'perexp_perfly',
      s = 'Mean mean fraction of frames per fly, per experiment';
  end
%   if donormbybinwidth,
%     s = [s,' / binwidth'];
%   end
  hylabel = ylabel(hax,s);
  
  % set ticks to reflect edge locations
  xtick = get(hax,'XTick');
  xticklabel = cellstr(get(hax,'XTickLabel'))';
  % this order makes sure that xtick will be chosen by unique instead of
  % edges
  xtick1 = [edges,xtick];
  xticklabel1 = [repmat({''},[1,nbins+1]),xticklabel];
  [xtick1,order] = unique(xtick1);
  xticklabel1 = xticklabel1(order);
  set(hax,'XTick',xtick1,'XTickLabel',xticklabel1,'TickDir','out',axesstyleparams{:});
  
  res.hxlabel = hxlabel;
  res.hylabel = hylabel;
  
  %% legend
  hlegend = legend(hax,legend_handles,legend_strings,'interpreter','none',legendstyleparams{:});
  res.hlegend = hlegend;
  
end



function my_frac = compute_frac(my_flies,my_ns)

  % take the intersection of the specified flies and experiments
  my_allflies_perexp = [movie2flies{my_ns}];
  my_flies = intersect(my_flies,my_allflies_perexp);
  my_ns = unique(fly2movie(my_flies));
  
  switch lower(averaging),
    case 'allexps_allflies',
      
      % treat all frames of data the same: just add up all the counts
      my_counts = sum(countsperfly(my_flies,:),1);
      my_frac = my_counts / sum(my_counts);
      
    case 'allexps_perfly',
      
      % get per-fly fracs, but treat all flies in all experiments the same
      my_frac = nanmean(fracperfly(my_flies,:),1);
      
    case 'perexp_allflies',
      
      % get per-exp fracs, but treat all frames within an experiment the
      % same
      
      % loop over experiments
      my_frac = zeros(1,nbins);
      my_nexps = zeros(1,nbins);
      for my_n = my_ns(:)',
        % which flies in this experiment
        my_flies_curr = intersect(my_flies,movie2flies{my_n});
        % frac for these flies
        my_frac_curr = sum(countsperfly(my_flies_curr,:),1);
        my_frac_curr = my_frac_curr / sum(my_frac_curr);
        % only include good data
        my_isdata = ~isnan(my_frac_curr);
        my_frac(my_isdata) = my_frac(my_isdata) + my_frac_curr(my_isdata);
        my_nexps(my_isdata) = my_nexps(my_isdata) + 1;
      end
      % normalize
      my_frac = my_frac ./ my_nexps;
      
    case 'perexp_perfly',
      
      % loop over experiments
      my_frac = zeros(1,nbins);
      my_nexps = zeros(1,nbins);
      for my_n = my_ns(:)',
        
        % which flies in this experiment
        my_flies_curr = intersect(my_flies,movie2flies{my_n});
        
        % average histogram for this experiment
        my_frac_curr = nanmean(fracperfly(my_flies_curr,:),1);
        
        % only include good data
        my_isdata = ~isnan(my_frac_curr);
        my_frac(my_isdata) = my_frac(my_isdata) + my_frac_curr(my_isdata);
        my_nexps(my_isdata) = my_nexps(my_isdata) + 1;
        
      end
      % normalize
      my_frac = my_frac ./ my_nexps;
      
    otherwise
      error('Unknown averaging method %s',averaging);
  end
  
end

end