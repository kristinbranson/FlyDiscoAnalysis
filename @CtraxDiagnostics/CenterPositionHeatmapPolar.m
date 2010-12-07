function heatmap = CenterPositionHeatmapPolar(obj,varargin)

%% parse inputs

% parse optional inputs
[edges_r,edges_theta,...
  lim_r,lim_theta,lim_prctile_theta,...
  nbins_r,nbins_theta,...
  conditions,...
  minspeed,...
  hax,hfig,doplot,...
  axesstyleparams,...
  ploterrorbars,plotperfly,plotperexp,clim,...
  figpos,fraclogscale,...
  squareaxes,...
  cm,...
  leftovers] = ...
  myparse_nocheck(varargin,...
  'edges_r',[],...
  'edges_theta',[],...
  'lim_r',nan(1,2),'lim_theta',nan(1,2),...
  'lim_prctile_theta',[1,99],...
  'nbins_r',obj.histogramtwomeasurements_nbins,...
  'nbins_theta',obj.histogramtwomeasurements_nbins,...
  'conditions',[],...
  'minspeed',0,...
  'hax',[],...
  'hfig',[],...
  'doplot',obj.histogramtwomeasurements_doplot,...
  'axesstyleparams',obj.histogramtwomeasurements_axesstyleparams,...
  'ploterrorbars',obj.histogramtwomeasurements_ploterrorbars,...
  'plotperfly',obj.histogramtwomeasurements_plotperfly,...
  'plotperexp',obj.histogramtwomeasurements_plotperexp,...
  'clim',nan(1,2),...
  'figpos',obj.histogramtwomeasurements_figpos,...
  'fraclogscale',false,...
  'squareaxes',true,...
  'colormap',jet(256));


%% set lim_r, lim_theta if not input

if isnan(lim_r(1)),
  lim_r(1) = 0;
end
if isnan(lim_r(2)),
  lim_r(2) = obj.arena_radius_mm;
end
if isnan(lim_theta(1)),
  lim_theta(1) = -pi;
end
if isnan(lim_theta(2)),
  lim_theta(2) = pi;
end

%% set edges in r if not input

% r edges: set binwidths so that the areas within each annulus is R^2 /
% nbins_r

% set edges if not input
if isempty(edges_r),
  r = lim_r(1); R = lim_r(2);
  edges_r = nan(1,nbins_r+1);
  edges_r(1) = r;
  for i = 2:nbins_r+1,
    dr = -r + sqrt(r^2 + R^2 / nbins_r);
    r = r + dr;
    edges_r(i) = r;
  end
end

%% add minspeed to conditions
if minspeed > 0,
  if isempty(conditions),
    conditions = @(trk) [false,trk.velmag > minspeed];
  else
    conditions = @(trk) conditions(trk) & [false,trk.velmag > minspeed];
  end
end

%% use 2-D histogram to compute, but not plot

heatmap = obj.HistogramTwoMeasurements('arena_r','wallangle',...
  'edges_x',edges_r,...
  'edges_y',edges_theta,'lim_y',lim_theta,'nbins_y',nbins_theta,'lim_prctile_y',lim_prctile_theta,...
  'conditions',conditions,...
  'doplot',false,leftovers{:});

ns = heatmap.ns;
flies = heatmap.flies;
nexpdirs = length(ns);
nflies = length(flies);
jackknife = heatmap.jackknife;
frac = heatmap.frac;
edges_theta = heatmap.edges_y;
averaging = heatmap.averaging;
expdirs = heatmap.expdirs;
movie2flies = heatmap.movie2flies;
if isfield(heatmap,'jackknife_frac_std'),
  jackknife_frac_std = heatmap.jackknife_frac_std;
end
if isfield(heatmap,'jackknife_frac_stderr'),
  jackknife_frac_stderr = heatmap.jackknife_frac_stderr;
end
fracperfly = heatmap.fracperfly;
fracperexp = heatmap.fracperexp;
Zperfly = heatmap.Zperfly;


%% plot

if doplot,
  
  %% create figures and axes

  % index within hfig
  nfigs = 1;
  figi_main = nfigs;
  doplotperexp = plotperexp && nexpdirs > 1;
  if doplotperexp,
    nfigs = nfigs + 1;
    figi_perexp = nfigs;
  end
  doplotperfly = plotperfly && nflies > 1;
  if doplotperfly,
    nfigs = nfigs + 1;
    figi_perfly = nfigs;
  end
  
  % index within hax
  nax = 1;
  axi_hist = nax;
  doploterrorbars = ~strcmpi(ploterrorbars,'none') && ~strcmpi(jackknife,'none');
  if doploterrorbars,
    nax = nax+1;
    axi_pluserror = nax;
    nax = nax+1;
    axi_minuserror = nax;
  end
  if doploterrorbars,
    nax_main = 3;
    axi_main = [axi_hist,axi_pluserror,axi_minuserror];
  else
    nax_main = 1;
    axi_main = axi_hist;
  end
  if doplotperexp,
    axi_perexp = nax + (1:nexpdirs);
    nax = nax + nexpdirs;
  end
  if doplotperfly,
    axi_perfly = nax + (1:nflies);
    nax = nax + nflies;
  end
  
  % get an axis for the main histogram, error hists
  % if hax does not exist for all of these, then create all of them
  if numel(hax) < max(axi_main) || ...
      ~all(ishandle(hax(axi_main))),
    % if no ax, check for figure
    if numel(hfig) < figi_main,
      % if no figure, create
      hfig(figi_main) = figure;
      if ~isempty(figpos),
        set(hfig(figi_main),'Position',figpos);
      end
    elseif ishandle(hfig(figi_main)),
      % there is a figure, so clear it
      clf(hfig(figi_main));
      figpos = get(hfig(figi_main),'Position');
    else
      % handle input, but figure does not exist, so make it
      figure(hfig(figi_main));
      if ~isempty(figpos),
        set(hfig(figi_main),'Position',figpos);
      end
    end
    % create the axes
    hax_main = createsubplots(1,nax_main,.05,hfig(figi_main));
    if doploterrorbars,
      hax_main = hax_main([2,1,3]);
    end
    hax(axi_main) = hax_main;
  else
    hfig(figi_main) = get(hax(axi_hist),'Parent');
  end
  
  if doplotperexp,
    % get axes for exps
    % if hax does not exist for all of these, then create all of them
    if ~isempty(figpos),
      figpos_perexp = get(hfig(figi_main),'Position');
      figpos_perexp(4) = figpos_perexp(4)*nax_main/nexpdirs;
    end
    if numel(hax) < max(axi_perexp) || ~all(ishandle(hax(axi_perexp))),
      % if no ax, check for figure
      if numel(hfig) < figi_perexp,
        % if no figure, create
        hfig(figi_perexp) = figure;
        if ~isempty(figpos),
          set(hfig(figi_perexp),'Position',figpos_perexp);
        end
      elseif ishandle(hfig(figi_perexp)),
        % there is a figure, so clear it
        clf(hfig(figi_perexp));
      else
        % handle input, but figure does not exist, so make it
        figure(hfig(figi_perexp));
        if ~isempty(figpos),
          set(hfig(figi_perexp),'Position',figpos_perexp);
        end
      end
      % create the axes
      hax(axi_perexp) = createsubplots(1,nexpdirs,.05,hfig(figi_perexp));
    else
      hfig(figi_perexp) = get(hax(axi_perexp(1)),'Parent');
    end
  end
  
  if doplotperfly,
    % get axes for flies
    % if hax does not exist for all of these, then create all of them
    if ~isempty(figpos),
      figpos_perfly = get(hfig(figi_main),'Position');
      figpos_perfly(4) = figpos_perfly(4)*nax_main/nflies;
    end
    if numel(hax) < max(axi_perfly) || ~all(ishandle(hax(axi_perfly))),
      % if no ax, check for figure
      if numel(hfig) < figi_perfly,
        % if no figure, create
        hfig(figi_perfly) = figure;
        if ~isempty(figpos),
          set(hfig(figi_perfly),'Position',figpos_perfly);
        end
      elseif ishandle(hfig(figi_perfly)),
        % there is a figure, so clear it
        clf(hfig(figi_perfly));
      else
        % handle input, but figure does not exist, so make it
        figure(hfig(figi_perfly));
        if ~isempty(figpos),
          set(hfig(figi_perfly),'Position',figpos_perfly);
        end
      end
      % create the axes
      hax(axi_perfly) = createsubplots(1,nflies,[[.01,.01];[.1,.1]],hfig(figi_perfly));
    else
      hfig(figi_perfly) = get(hax(axi_perfly(1)),'Parent');
    end
  end
  
  % clear the axes
  for axi = 1:nax,
    cla(hax(axi));
  end  

  % store
  heatmap.plotparams.hfig = hfig;
  heatmap.plotparams.hax = hax;
  
  % for clim
  minv = inf;
  maxv = -inf;
  
  %% plot the main histogram
  hhist = polarimagesc(edges_r,edges_theta,obj.arena_center_mm,frac',...
    'parent',hax(axi_hist),'tag','mainhist');
  colorbar('peer',hax(axi_hist));
  % update limits
  minv = min(minv,min(frac(~isinf(frac))));
  maxv = max(maxv,max(frac(~isinf(frac))));

  heatmap.plotparams.hhist = hhist;

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
  
  title(hax(axi_hist),s);

  % user data for main plot
  plotinfo = struct;
  plotinfo.expdirs = expdirs;
  plotinfo.flies = flies;
  plotinfo.averaging = averaging;
  plotinfo.jackknife = jackknife;
  set(hax(axi_hist),'UserData',plotinfo);
  
  %% plot standard error/standard deviation
  
  if doploterrorbars,
    switch ploterrorbars,
      case 'std',
        err = jackknife_frac_std;
      case 'stderr',
        err = jackknife_frac_stderr;
      otherwise
        error('Unknown ploterrorbars type %s',ploterrorbars);
    end
    
    errinfo = plotinfo;
    errinfo.jackknife = jackknife;
    if strcmpi(jackknife,'perfly'),
      errinfo.n = nflies;
    elseif strcmpi(jackknife,'perexp'),
      errinfo.n = nexpdirs;
    end
    
    im = max(0,frac-err);
    herr_minus = polarimagesc(edges_r,edges_theta,obj.arena_center_mm,im',...
      'parent',hax(axi_minuserror),'tag',sprintf('minuserr_%s',ploterrorbars),...
      'userdata',errinfo);
    colorbar('peer',hax(axi_minuserror));
    % update limits
    minv = min(minv,min(im(~isinf(im))));

    im = frac+err;
    herr_plus = polarimagesc(edges_r,edges_theta,obj.arena_center_mm,im',...
      'parent',hax(axi_pluserror),'tag',sprintf('pluserr_%s',ploterrorbars),...
      'userdata',errinfo);
    colorbar('peer',hax(axi_pluserror));
    
    maxv = max(maxv,max(im(~isinf(im))));

    heatmap.plotparams.herr_minus = herr_minus;
    heatmap.plotparams.herr_plus = herr_plus;
    
    s = sprintf('- %s',ploterrorbars);
    title(hax(axi_minuserror),s);
    s = sprintf('+ %s',ploterrorbars);
    title(hax(axi_pluserror),s);
        
  end
  
  %% plot per-experiment histogram
  
  if doplotperexp,
    hperexphist = zeros(1,nexpdirs);
    for n = 1:nexpdirs,
      expinfo = struct('n',ns(n),'flies',flies(movie2flies{n}),'expdir',expdirs{n},'ndataperfly',Zperfly(movie2flies{n}));
      fracperexpcurr = permute(fracperexp(n,:,:),[2,3,1]);
      im = fracperexpcurr;

      % update limits
      minv = min(minv,min(im(~isinf(im))));
      maxv = max(maxv,max(im(~isinf(im))));
      hperexphist(n) = polarimagesc(edges_r,edges_theta,obj.arena_center_mm,im',...
        'parent',hax(axi_perexp(n)),'tag',sprintf('exp%d',n),'userdata',expinfo);
      title(hax(axi_perexp(n)),expdirs{n},'interpreter','none');
    end
    heatmap.plotparams.hperexphist = hperexphist;
    
  end
  
  %% plot the per-fly histograms
  if plotperfly && nflies > 1,
    hperflyhist = zeros(1,nflies);
    for fly = 1:nflies,
      flyglobal = flies(fly);
      nglobal = obj.fly2movie(flyglobal);
      expdir = obj.expdir_bases{nglobal};
      flyofmovie = find(flyglobal == obj.movie2flies{nglobal},1);
      flyinfo = struct('fly',flyglobal,'n',nglobal,'expdir',expdir,'ndata',Zperfly(fly),'flyofmovie',flyofmovie);
      
      fracperflycurr = permute(fracperfly(fly,:,:),[2,3,1]);
      im = fracperflycurr;

      hperflyhist(fly) = polarimagesc(edges_r,edges_theta,obj.arena_center_mm,im',...
        'parent',hax(axi_perfly(fly)),'tag',sprintf('fly%d',fly),'userdata',flyinfo);
      
      title(hax(axi_perfly(fly)),sprintf('%d,%d',nglobal,flyglobal),'interpreter','none');
      axis(hax(axi_perfly(fly)),'off');
    end
    
    % update limits
    minv = min(minv,min(im(~isinf(im))));
    maxv = max(maxv,max(im(~isinf(im))));

    heatmap.plotparams.hperflyhist = hperflyhist;
    
  end
  
  %% set axis limits

  if isnan(clim(1)),
    clim(1) = minv;
  end
  if isnan(clim(2)),
    clim(2) = maxv;
  end
  for axi = 1:numel(hax),
    set(hax(axi),'clim',clim);
    axis(hax(axi),'xy','tight');
    if squareaxes,
      axis(hax(axi),'square');
    end
    if ~isempty(axesstyleparams),
      set(hax(axi),axesstyleparams{:});
    end
  end
  linkaxes(hax(1:nax));
  linkprop(hax(1:nax),'clim');
  
  heatmap.plotparams.clim = clim;
  
  %% set colormap

  if fraclogscale,
    cm = logscale_colormap(cm,clim);
  end
  for figi = 1:nfigs
    set(hfig(figi),'Colormap',cm);
  end
  
  %% labels
  
  % title
  plottitle = 'Center position heatmap, polar';
  set(hfig(figi_main),'name',plottitle);
  if doplotperexp,
    set(hfig(figi_perexp),'name',[plottitle,' per exp']);
  end
  if doplotperfly,
    set(hfig(figi_perfly),'name',[plottitle,' per fly']);
  end
      
end