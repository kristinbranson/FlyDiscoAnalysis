function res = HistogramTwoMeasurements(obj,fn_x,fn_y,varargin)

res = struct;

%% parse inputs

% check that there is at least one experiment loaded
if obj.nexpdirs == 0,
  error('No experiments loaded');
end

% check that the input measurements are defined
if ~isfield(obj.trx,fn_x),
  error('Field %s not computed for trx',fn_x);
end
if ~isfield(obj.trx,fn_y),
  error('Field %s not computed for trx',fn_y);
end

% parse optional inputs
[edges_x,edges_y,...
  nbins_x,nbins_y,nbins,...
  lim_x,lim_y,...
  lim_prctile_x,lim_prctile_y,...
  hax,hfig,doplot,...
  conditions,...
  outputfun_x,outputfun_string_x,outputfun_y,outputfun_string_y,...
  expdirs,flies,averaging,jackknife,...
  axesstyleparams,...
  ploterrorbars,plotperfly,plotperexp,clim,docomputestd,docomputestderr,...
  binmode_x,binmode_y,figpos,fraclogscale,...
  squareaxes,...
  cm] = ...
  myparse(varargin,...
  'edges_x',[],'edges_y',[],...
  'nbins_x',obj.histogramtwomeasurements_nbins_x,...
  'nbins_y',obj.histogramtwomeasurements_nbins_y,...
  'nbins',obj.histogramtwomeasurements_nbins,...
  'lim_x',nan(1,2),'lim_y',nan(1,2),...
  'lim_prctile_x',[1,99],'lim_prctile_y',[1,99],...
  'hax',[],...
  'hfig',[],...
  'doplot',obj.histogramtwomeasurements_doplot,...
  'conditions',[],...
  'outputfun_x',[],...
  'outputfun_string_x','',...
  'outputfun_y',[],...
  'outputfun_string_y','',...
  'expdirs',obj.expdir_bases,...
  'flies',1:obj.nflies,...
  'averaging',obj.histogramtwomeasurements_averaging,...
  'jackknife',obj.histogramtwomeasurements_jackknife,...
  'axesstyleparams',obj.histogramtwomeasurements_axesstyleparams,...
  'ploterrorbars',obj.histogramtwomeasurements_ploterrorbars,...
  'plotperfly',obj.histogramtwomeasurements_plotperfly,...
  'plotperexp',obj.histogramtwomeasurements_plotperexp,...
  'clim',nan(1,2),...
  'docomputestd',obj.histogramtwomeasurements_docomputestd,...
  'docomputestderr',obj.histogramtwomeasurements_docomputestderr,...
  'binmode_x','linear',...
  'binmode_y','linear',...
  'figpos',obj.histogramtwomeasurements_figpos,...
  'fraclogscale',false,...
  'squareaxes',true,...
  'colormap',jet(256));

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

% if jackknife not set, then choose perexp if there is more than one
% experiment, else choose perfly if there is more than one fly, else choose
% none
jackknife = SetJackKnifeMethod(jackknife,nflies,nexpdirs);
% store for returning
res.jackknife = jackknife;
  
%% set edges if not input

% x edges
% nbins_x defaults to nbins
% nbins should not be nan
if isempty(nbins_x),
  nbins_x = nbins;
end
[edges_x,nbins_x,centers_x] = obj.SelectHistEdges(fn_x,edges_x,flies,conditions,nbins_x,lim_x,lim_prctile_x,outputfun_x,binmode_x);

% y edges
% nbins_y defaults to nbins
% nbins should not be nan
if isempty(nbins_y),
  nbins_y = nbins;
end
[edges_y,nbins_y,centers_y] = obj.SelectHistEdges(fn_y,edges_y,flies,conditions,nbins_y,lim_y,lim_prctile_y,outputfun_y,binmode_y);

% store for returning
res.edges_x = edges_x;
res.centers_x = centers_x;
res.nbins_x = nbins_x;
res.edges_y = edges_y;
res.centers_y = centers_y;
res.nbins_y = nbins_y;

%% histogram for each fly

countsperfly = nan([length(flies),nbins_y,nbins_x]);
for i = 1:nflies,
  
  % current fly
  fly = flies(i);
  
  % get x and y data for this fly
  if ~isempty(conditions),
    
    % filter frames by conditions
    idx = conditions(obj.trx(fly));
    x = obj.trx(fly).(fn_x)(idx);
    y = obj.trx(fly).(fn_y)(idx);

  else
    x = obj.trx(fly).(fn_x);
    y = obj.trx(fly).(fn_y);
  end
  
  % apply outputfun to x, y
  if ~isempty(outputfun_x),
    x = outputfun_x(x);
  end
  if ~isempty(outputfun_y),
    y = outputfun_y(y);
  end
  
  % actual histogramming!
  countscurr = hist3([y(:),x(:)],'Edges',{edges_y,edges_x});
  % last value is the number that equal the right-most edge, combine with
  % last bin
  countscurr(:,end-1) = countscurr(:,end-1)+countscurr(:,end);
  countscurr(end-1,:) = countscurr(end-1,:)+countscurr(end,:);
  countscurr = countscurr(1:end-1,1:end-1);
  % store
  countsperfly(i,:,:) = countscurr;
end

%% set up to index into data structures only computed for selected flies, selected exps

[expdirs,movie2flies,fly2movie] = ...
  obj.SubsetDataStructs(flies,ns);

res.expdirs = expdirs;
res.movie2flies = movie2flies;
res.fly2movie = fly2movie;

%% per-fly histograms

% these won't all be used for all averaging conditions, but it should be
% quick and will speed up jackknifing, so compute them
Zperfly = sum(sum(countsperfly,2),3);
fracperfly = bsxfun(@rdivide,countsperfly,Zperfly);
res.fracperfly = fracperfly;
res.Zperfly = Zperfly;

%% per-exp histograms

fracperexp = nan(nexpdirs,nbins_y,nbins_x);
for n = 1:nexpdirs,
  fracperexp(n,:,:) = CtraxStatsBase.CollateHistograms(1:nflies,n,...
    countsperfly,movie2flies,fly2movie,averaging,fracperfly);
end
res.fracperexp = fracperexp;

%% compute frac using all the data

% compute final histogram using all the data
frac = CtraxStatsBase.CollateHistograms(1:nflies,1:nexpdirs,...
  countsperfly,movie2flies,fly2movie,averaging,fracperfly);

res.frac = frac;
res.averaging = averaging;

%% jackknife to get an estimate of standard deviation

jackknife_fun_handle = @(flies_in,ns_in) CtraxStatsBase.CollateHistograms(flies_in,ns_in,...
  countsperfly,movie2flies,fly2movie,averaging,fracperfly);

if docomputestd && ~strcmpi(jackknife,'none'),
  
  jackknife_frac_std = ...
    CtraxStatsBase.JackKnifeStd(jackknife_fun_handle,nflies,nexpdirs,jackknife);

end

% store for returning
res.jackknife_frac_std = jackknife_frac_std;

%% jackknife to get an estimate of standard error
if docomputestderr && ~strcmpi(jackknife,'none'),

  jackknife_frac_stderr = ...
    CtraxStatsBase.JackKnifeStdErr(jackknife_fun_handle,nflies,nexpdirs,jackknife);
   
end

% store for returning
res.jackknife_frac_stderr = jackknife_frac_stderr;
    
%% plot

res.plotparams = struct;

if doplot,

  doplotperexp = plotperexp && nexpdirs > 1;
  doplotperfly = plotperfly && nflies > 1;
  doploterrorbars = ~strcmpi(ploterrorbars,'none') && ~strcmpi(jackknife,'none');
  docla = true;
  
  %% create figures and axes
  [hfig,hax,nfigs,figi_main,figi_perexp,figi_perfly,...
    nax,axi_hist,axi_pluserror,axi_minuserror,axi_perexp,axi_perfly] = ...
    CtraxStatsBase.Create2DHistogramFigures(hax,hfig,nflies,nexpdirs,...
    doplotperexp,doplotperfly,doploterrorbars,figpos,docla);
  
  % store
  res.plotparams.hfig = hfig;
  res.plotparams.hax = hax;
  
  % for clim
  minv = inf;
  maxv = -inf;
  
  %% plot the main histogram
  im = frac;
  hhist = imagesc(edges_x([1,nbins_x+1]),...
    edges_y([1,nbins_y+1]),im,'parent',hax(axi_hist),...
    'tag','mainhist');
  colorbar('peer',hax(axi_hist));
  % update limits
  minv = min(minv,min(im(~isinf(im))));
  maxv = max(maxv,max(im(~isinf(im))));

  res.plotparams.hhist = hhist;

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
    herr_minus = imagesc(edges_x([1,nbins_x+1]),...
      edges_y([1,nbins_y+1]),im,'parent',hax(axi_minuserror),...
      'tag',sprintf('minuserr_%s',ploterrorbars),'userdata',errinfo);
    colorbar('peer',hax(axi_minuserror));
    % update limits
    minv = min(minv,min(im(~isinf(im))));

    im = frac+err;
    herr_plus = imagesc(edges_x([1,nbins_x+1]),...
      edges_y([1,nbins_y+1]),im,'parent',hax(axi_pluserror),...
      'tag',sprintf('pluserr_%s',ploterrorbars),'userdata',errinfo);
    colorbar('peer',hax(axi_pluserror));
    
    maxv = max(maxv,max(im(~isinf(im))));

    res.plotparams.herr_minus = herr_minus;
    res.plotparams.herr_plus = herr_plus;
    
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
      hperexphist(n) = ...
        imagesc(edges_x([1,nbins_x+1]),...
        edges_y([1,nbins_y+1]),im,...
        'parent',hax(axi_perexp(n)),...
        'tag',sprintf('exp%d',n),'userdata',expinfo);

      title(hax(axi_perexp(n)),expdirs{n}(end-20:end),'interpreter','none');
    end
    res.plotparams.hperexphist = hperexphist;
    
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

      hperflyhist(fly) = ...
        imagesc(edges_x([1,nbins_x+1]),...
        edges_y([1,nbins_y+1]),im,...
        'parent',hax(axi_perfly(fly)),...
        'tag',sprintf('fly%d',fly),'userdata',flyinfo);
      title(hax(axi_perfly(fly)),sprintf('%d,%d',nglobal,flyglobal),'interpreter','none');
      axis(hax(axi_perfly(fly)),'off');
    end
    
    % update limits
    minv = min(minv,min(im(~isinf(im))));
    maxv = max(maxv,max(im(~isinf(im))));

    res.plotparams.hperflyhist = hperflyhist;
    
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
    axis(hax(axi),'xy');
    if squareaxes,
      axis(hax(axi),'square');
    end
    if ~isempty(axesstyleparams),
      set(hax(axi),axesstyleparams{:});
    end
  end
  linkaxes(hax(1:nax));
  linkprop(hax(1:nax),'clim');
  
  res.plotparams.clim = clim;
  
  %% set colormap

  if fraclogscale,
    cm = logscale_colormap(cm,clim);
  end
  for figi = 1:nfigs
    set(hfig(figi),'Colormap',cm);
  end

  %% labels

  % x, y labels
  if ~isempty(outputfun_x),
    if ~isempty(outputfun_string_x),
      s_x = sprintf('%s of %s',outputfun_string_x,fn_x);
    else
      s_x = sprintf('%s of %s',func2str(outputfun_x),fn_x);
    end
  else
    s_x = fn_x;
  end
  if ~isempty(outputfun_y),
    if ~isempty(outputfun_string_y),
      s_y = sprintf('%s of %s',outputfun_string_y,fn_y);
    else
      s_y = sprintf('%s of %s',func2str(outputfun_y),fn_y);
    end
  else
    s_y = fn_y;
  end
  hxlabel = nan(size(hax));
  hylabel = nan(size(hax));
  for axi = 1:nax,
    hxlabel(axi) = xlabel(hax(axi),s_x,'interpreter','none');
    hylabel(axi) = ylabel(hax(axi),s_y,'interpreter','none');
  end
  
  % title
  plottitle = sprintf('%s vs %s',s_x,s_y);
  set(hfig(figi_main),'name',plottitle);
  if doplotperexp,
    set(hfig(figi_perexp),'name',[plottitle,' per exp']);
  end
  if doplotperfly,
    set(hfig(figi_perfly),'name',[plottitle,' per fly']);
  end
  
  res.hxlabel = hxlabel;
  res.hylabel = hylabel;
    
end
