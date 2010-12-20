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
  
  doplotperexp = plotperexp && nexpdirs > 1;
  doplotperfly = plotperfly && nflies > 1;
  doploterrorbars = ~strcmpi(ploterrorbars,'none') && ~strcmpi(jackknife,'none');
  docla = true;
  
  %% create figures and axes
  [hfig,hax,nfigs,figi_main,figi_perexp,figi_perfly,...
    nax,axi_hist,axi_pluserror,axi_minuserror,axi_perexp,axi_perfly] = ...
    CraxStatsBase.Create2DHistogramFigures(hax,hfig,nflies,nexpdirs,...
    doplotperexp,doplotperfly,doploterrorbars,figpos,docla);
  
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
      title(hax(axi_perexp(n)),expdirs{n}(end-20:end),'interpreter','none');
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