function res = BowlAngleBias(obj,varargin)

%% parse inputs

[nbins,...
  anglerange,...
  min_radius,...
  hfig,...
  hax,...
  doplot,...
  minspeed,...
  conditions,...
  expdirs,...
  flies,...
  averaging,...
  jackknife,...
  axesstyleparams,...
  errstyleparams,...
  perflyhiststyleparams,...
  perexphiststyleparams,...
  legendstyleparams,...
  ploterrorbars,plotperfly,plotperexp,...
  docomputestd,docomputestderr,...
  figpos,...
  axisplot,...
  trange] = ...
  myparse(varargin,...
  'nbins',360,...
  'anglerange',pi/2,...
  'min_radius',0,...
  'hfig',[],...
  'hax',[],...
  'doplot',obj.histogrammeasurements_doplot,...
  'minspeed',0,...
  'conditions',[],...
  'expdirs',obj.expdir_bases,...
  'flies',1:obj.nflies,...
  'averaging',obj.histogrammeasurements_averaging,...
  'jackknife',obj.histogrammeasurements_jackknife,...
  'axesstyleparams',obj.histogrammeasurements_axesstyleparams,...
  'errstyleparams',{},...
  'perflyhiststyleparams',obj.histogrammeasurements_perflyhiststyleparams,...
  'perexphiststyleparams',obj.histogrammeasurements_perexphiststyleparams,...
  'legendstyleparams',obj.histogrammeasurements_legendstyleparams,...
  'ploterrorbars',obj.histogrammeasurements_ploterrorbars,...
  'plotperfly',obj.histogrammeasurements_plotperfly,...
  'plotperexp',obj.histogrammeasurements_plotperexp,...
  'docomputestd',obj.histogrammeasurements_docomputestd,...
  'docomputestderr',obj.histogrammeasurements_docomputestderr,...
  'figpos',obj.histogrammeasurements_figpos,...
  'axisplot',[],...
  'trange',nan(1,2));

res.anglerange = anglerange;
res.min_radius = min_radius;
res.minspeed = minspeed;
res.trange = trange;

%% add minspeed to conditions
if minspeed > 0 || min_radius > 0,
  if isempty(conditions),
    conditions = @(trk) [false,(trk.velmag >= minspeed)] & (trk.arena_r >= min_radius);
  else
    conditions = @(trk) conditions(trk) & [false,trk.velmag >= minspeed] & ...
      (trk.arena_r >= min_radius);
  end
end

%% intersect flies and experiments

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

%% edges
edges = linspace(-pi,pi,nbins+1);
bin_width = edges(2)-edges(1);
bin_centers = (edges(1:end-1)+edges(2:end))/2;

%% normalize so that the areas in the two sides are the same

% how many bins fall within range
nbinsin = floor(anglerange / bin_width);
fil = ones(1,nbinsin);

% center of the range
% nin(nbinsin) = sum(counts(1:nbinsin))
% so thetares(nbinsin) = (theta(1)+theta(nbinsin))/2
% = (theta(nbinsin) - bin_width*(nbinsin-1) + theta(nbinsin)) /2
% = theta(nbinsin) - bin_width*(nbinsin-1)/2
theta_centers = modrange(bin_centers - bin_width*(nbinsin-1)/2,-pi,pi);

% normalize so that the areas in the two sides are the same
fracbinsin = nbinsin / nbins;

res.fracbinsin = fracbinsin;

%% histogram data into many bins

countsperfly = nan(length(flies),nbins);
for i = 1:nflies,
  fly = flies(i);
  if ~isempty(conditions),
    x = obj.trx(fly).wallangle(conditions(obj.trx(fly)));
  else
    x = obj.trx(fly).wallangle;
  end
  countscurr = histc(x,edges);
  countsperfly(i,:) = [countscurr(1:end-2),countscurr(end-1)+countscurr(end)];
end

Zperfly = sum(countsperfly,2);

res.countsperfly = countsperfly;
res.Zperfly = Zperfly;

%% set up to index into data structures only computed for selected flies, selected exps

[expdirs,movie2flies,fly2movie,nfliespermovie] = ...
  obj.SubsetDataStructs(flies,ns);

%% compute bias using all data

bias_fun_handle = @(f) compute_bias(f,nbins,fracbinsin,fil);

[bias,biasperfly] = ...
  CtraxDiagnostics.CollateHistStats(bias_fun_handle,...
  1:nflies,1:nexpdirs,countsperfly,movie2flies,fly2movie,averaging);

res.bias = bias;
res.biasperfly = biasperfly;

%% jackknife to get an estimate of standard deviation

jackknife_fun_handle = @(flies_in,ns_in) ...
  CtraxDiagnostics.CollateHistStats(bias_fun_handle,flies_in,ns_in,...
  countsperfly,movie2flies,fly2movie,averaging,biasperfly);

if docomputestd && ~strcmpi(jackknife,'none'),
  
  jackknife_bias_std = ...
    CtraxDiagnostics.JackKnifeStd(jackknife_fun_handle,nflies,nexpdirs,jackknife);
  res.jackknife_bias_std = jackknife_bias_std;

end

%% jackknife to get an estimate of standard error

if docomputestderr && ~strcmpi(jackknife,'none'),
  
  jackknife_bias_stderr = ...
    CtraxDiagnostics.JackKnifeStdErr(jackknife_fun_handle,nflies,nexpdirs,jackknife);
  res.jackknife_bias_stderr = jackknife_bias_stderr;

end

%% plot

if doplot
  
  %% get axes to plot in
  [hax,hfig] = get_axes(hax,hfig,'figpos',figpos);
  cla(hax);
  hold(hax,'on');
  
  res.plotparams = struct;
  
  res.plotparams.hax = hax;
  res.plotparams.hfig = hfig;
  
  %% colors for positive, negative bias, error bar region
  poscolor = [.7,0,0];
  negcolor = [0,0,1];
  errlighten = .7;
  
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
  
  % plot positive and negative biases separatly
  ispos = bias > 0;
  isneg = bias < 0;
  
  if ~strcmpi(ploterrorbars,'none') && ~strcmpi(jackknife,'none')
    switch ploterrorbars,
      case 'std',
        err = jackknife_bias_std;
      case 'stderr',
        err = jackknife_bias_stderr;
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
    
    % what color is the error region
    %ispos_strict = bias - err >= 0;
    %isneg_strict = bias + err <= 0;
    %isneither_strict = ~ispos_strict & ~isneg_strict;
    
    herr = nan(1,4);
    
    % use bins to get continuity
    theta_err = [theta_centers+bin_width/2,theta_centers(1)-bin_width/2];
    bias_err = [bias,bias(1)];
    err_err = [err,err(1)];
    ispos_err = [false,ispos] | [ispos,false];
    isneg_err = [false,isneg] | [isneg,false];
    
    % positive
    r_plot = abs(bias_err);
    r_plot(~ispos_err) = 0;
    r_err = err_err;
    r_err(~ispos_err) = 0;
    herr(1) = polarerrorpatch(theta_err,r_plot,r_err,poscolor*(1-errlighten)+errlighten,...
      'parent',hax,'edgecolor','none','tag',sprintf('pos_err_%s',ploterrorbars),...
      'userdata',errinfo);
    
    % negative
    r_plot = abs(bias_err);
    r_plot(~isneg_err) = 0;
    r_err = err_err;
    r_err(~isneg_err) = 0;
    herr(2) = polarerrorpatch(theta_err,r_plot,r_err,negcolor*(1-errlighten)+errlighten,...
      'parent',hax,'edgecolor','none','tag',sprintf('neg_err_%s',ploterrorbars),...
      'userdata',errinfo);
    
    % negative error bound, positive bias
    idx = ispos_err & (bias_err - err_err < 0);
    r_plot = zeros(size(bias_err));
    r_err = err_err - bias_err;
    r_err(~idx) = 0;
    herr(3) = polarerrorpatch(theta_err,r_plot,r_err,(poscolor+negcolor)/2*(1-errlighten)+errlighten,...
      'parent',hax,'edgecolor','none','tag',sprintf('pos_err_overlap_%s',ploterrorbars),...
      'userdata',errinfo);
    
    % positive error bound, negative bias
    idx = isneg_err & (bias_err + err_err > 0);
    r_plot = zeros(size(bias_err));
    r_err = err_err + bias_err;
    r_err(~idx) = 0;
    herr(4) = polarerrorpatch(theta_err,r_plot,r_err,(poscolor+negcolor)/2*(1-errlighten)+errlighten,...
      'parent',hax,'edgecolor','none','tag',sprintf('neg_err_overlap_%s',ploterrorbars),...
      'userdata',errinfo);
    if ~isempty(errstyleparams),
      set(herr,errstyleparams{:});
    end
    
    legend_handles(end+1) = herr(3);
    legend_strings{end+1} = sprintf('%s (bagging %s)',ploterrorbars,jackknife);

    res.plotparams.herr = herr;
  
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
  
  %% plot the per-fly histograms
  if plotperfly,
    hperflyhist = zeros(1,nflies);
    for fly = 1:nflies,
      flyglobal = flies(fly);
      nglobal = obj.fly2movie(flyglobal);
      expdir = obj.expdir_bases{nglobal};
      flyofmovie = find(flyglobal == obj.movie2flies{nglobal},1);
      flyinfo = struct('fly',flyglobal,'n',nglobal,'expdir',expdir,'ndata',Zperfly(fly),'flyofmovie',flyofmovie);

      x_plot = abs(biasperfly(fly,:)) .* cos(theta_centers);
      y_plot = abs(biasperfly(fly,:)) .* sin(theta_centers);
      
      hperflyhist(fly) = plot(hax,x_plot,y_plot,'color',colorsperfly(fly,:),...
        'linewidth',.5,'userdata',flyinfo,'tag',sprintf('perfly%d',fly));
    end
    
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
      biasperexpcurr = ...
        CtraxDiagnostics.CollateHistStats(bias_fun_handle,1:nflies,n,...
        countsperfly,movie2flies,fly2movie,averaging,biasperfly);
      x_plot = abs(biasperexpcurr) .* cos(theta_centers);
      y_plot = abs(biasperexpcurr) .* sin(theta_centers);
      hperexphist(n) = plot(hax,x_plot,y_plot,'-','color',colorsperexp(n,:),...
        'linewidth',1,'tag',sprintf('exp%d',n),'userdata',expinfo);
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
    
  
  %% plot the average bias
  
  x_plot = abs(bias) .* cos(theta_centers);
  y_plot = abs(bias) .* sin(theta_centers);
  
  hbias = nan(1,2);
  
  x_plot_curr = x_plot; y_plot_curr = y_plot;
  x_plot_curr(~ispos) = nan;
  y_plot_curr(~ispos) = nan;
  hbias(1) = plot(x_plot_curr,y_plot_curr,'-','color',poscolor,'linewidth',2);
  x_plot_curr = x_plot; y_plot_curr = y_plot;
  x_plot_curr(~isneg) = nan;
  y_plot_curr(~isneg) = nan;
  hbias(2) = plot(x_plot_curr,y_plot_curr,'-','color',negcolor,'linewidth',2);
  
  res.plotparams.hhist = hbias;
  
  legend_handles(end+1) = hbias(1);
  legend_strings{end+1} = 'Positive bias, all flies, exps';
  legend_handles(end+1) = hbias(2);
  legend_strings{end+1} = 'Negative bias, all flies, exps';
  
  %% set axis limits
  if isempty(axisplot),
    axisplot = [-1,1,-1,1]/fracbinsin;
  end
  axis(hax,'equal');
  axis(hax,axisplot);
  
  res.plotparams.axisplot = axisplot;
  
  if ~isempty(axesstyleparams),
    set(hax,axesstyleparams{:});
  end

  
  %% labels
  
  switch averaging,
    case 'allexps_allflies',
      s = 'Bias over all frames';
    case 'allexps_perfly',
      s = 'Mean bias per fly';
    case 'perexp_allflies',
      s = 'Mean bias over all frames per experiment';
    case 'perexp_perfly',
      s = 'Mean mean bias per fly, per experiment';
  end
  s = sprintf('%s, min speed = %.1f',s,minspeed);
  if all(~isnan(trange)),
    s = sprintf('%s, t in [%.1f,%.1f]',s,trange(1),trange(2));
  elseif ~isnan(trange(1)),
    s = sprintf('%s, t >= %.1f',s,trange(1));
  elseif ~isnan(trange(2)),
    s = sprintf('%s, t <= %.1f',s,trange(2));
  end
  
  title(hax,s);
  set(hfig,'name','Bias');
    
  %% legend
  hlegend = legend(hax,legend_handles,legend_strings,'interpreter','none',legendstyleparams{:});
  res.plotparams.hlegend = hlegend;
  
end

function bias = compute_bias(fracperbin,nbins,fracbinsin,fil)

% bias = fracin / fracbinsin - fracout / fracbinsout;

% circular filter to get sum in each sequence of nbinsin bins
% nin(nbinsin) = sum(counts(1:nbinsin))
bias = nan(size(fracperbin));
for i = 1:size(fracperbin,1),
  fracin = cconv(fil,fracperbin(i,:),nbins);
  fracout = 1 - fracin;
  bias(i,:) = fracin / fracbinsin - fracout / (1 - fracbinsin);
end
