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

%% take the intersection of specified flies and expdirs

[ns,flies] = obj.IntersectFliesExpdirs(flies,expdirs);
nflies = length(flies);
nexpdirs = length(ns);

% save for returning
res.ns = ns;
res.flies = flies;

% check that there is at least one experiment selected
if nexpdirs == 0,
  error('No experiments selected.');
end

%% set jackknife method if not set

jackknife = SetJackKnifeMethod(jackknife,nflies,nexpdirs);
res.jackknife = jackknife;
  
%% set edges if not input

[edges,nbins,centers] = SelectHistEdges(obj,...
  fn,edges,flies,conditions,nbins,lim,lim_prctile,outputfun,binmode);

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

[expdirs,movie2flies,fly2movie,nfliespermovie] = ...
  obj.SubsetDataStructs(flies,ns);

%% per-fly histograms

% these won't all be used for all averaging conditions, but it should be
% quick and will speed up jackknifing, so compute them
Zperfly = sum(countsperfly,2);
fracperfly = bsxfun(@rdivide,countsperfly,Zperfly);

%% compute frac using all the data

% compute histogram using all the data
frac = CtraxStatsBase.CollateHistograms(1:nflies,1:nexpdirs,...
  countsperfly,movie2flies,...
  fly2movie,averaging,fracperfly);

res.frac = frac;

%% jackknife to get an estimate of standard deviation

jackknife_fun_handle = @(flies_in,ns_in) CtraxStatsBase.CollateHistograms(flies_in,ns_in,...
  countsperfly,movie2flies,fly2movie,averaging,fracperfly);

if docomputestd && ~strcmpi(jackknife,'none'),

  jackknife_frac_std = ...
    CtraxStatsBase.JackKnifeStd(jackknife_fun_handle,nflies,nexpdirs,jackknife);
  
end

res.jackknife_frac_std = jackknife_frac_std;

%% jackknife to get an estimate of standard error

if docomputestderr && ~strcmpi(jackknife,'none'),

  jackknife_frac_stderr = ...
    CtraxStatsBase.JackKnifeStdErr(jackknife_fun_handle,nflies,nexpdirs,jackknife);
  
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
  [hax,hfig] = get_axes(hax,hfig,'figpos',figpos);
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
    legendexps = 1:nexpdirs;
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
    
    if ylogscale,
      minv = min(frac(:))/2;
    else
      minv = 0;
    end
    herr = patch([centers,fliplr(centers)],max(minv,[frac-err,fliplr(frac+err)]),.8*ones(1,3),...
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
      fracperexpcurr = CtraxStatsBase.CollateHistograms(1:nflies,n,...
        countsperfly,movie2flies,fly2movie,averaging,fracperfly);
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
