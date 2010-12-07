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
  squareaxes] = ...
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
  'squareaxes',true);

%% take the intersection of specified flies and expdirs

% which expdirs will we look at?
[didfind,ns] = ismember(expdirs,obj.expdir_bases);
% report if requested expdirs are not loaded
if ~all(didfind),
  warning(['The following expdirs are not loaded: ',sprintf('%s ',expdirs{~didfind})]);
end
ns = ns(didfind);

% which flies does this correspond to?
allflies_perexp = [obj.movie2flies{ns}];
% take the intersection of the specified flies and the flies in the
% specified experiments
flies = intersect(flies,allflies_perexp);
% which experiments do these flies correspond to
ns = unique(obj.fly2movie(flies));
nflies = length(flies);
nexpdirs = length(ns);

% save for returning
res.ns = ns;
res.nflies = nflies;

% check that there is at least one experiment selected
if nexpdirs == 0,
  error('No experiments selected.');
end

%% set jackknife method if not set

% if jackknife not set, then choose perexp if there is more than one
% experiment, else choose perfly if there is more than one fly, else choose
% none
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
% if there isn't set to perfly
if nexpdirs == 1 && strcmpi(jackknife,'perexp'),
  warning('Only one experiment selected, but jackknife = ''perexp''. Splitting over flies instead.');
  jackknife = 'perfly';
end
% to split over flies, we need more than one fly. if there isn't set to
% none
if nflies == 1 && strcmpi(jackknife,'perfly'),
  warning('Only one fly selected, but jackknife = ''perfly''. Splitting over frames not implemented. Not jackknifing.');
  jackknife = 'none';
end

% store for returning
res.jackknife = jackknife;
  
%% set edges if not input

% x edges
% nbins_x defaults to nbins
% nbins should not be nan
if isempty(nbins_x),
  nbins_x = nbins;
end
% we'll need to do this twice, so put it in a function:
% choose edges_x
edges_x = get_edges(fn_x,edges_x,nbins_x,lim_x,lim_prctile_x,outputfun_x,binmode_x);
centers_x = (edges_x(1:end-1)+edges_x(2:end))/2;

% y edges
% nbins_y defaults to nbins
% nbins should not be nan
if isempty(nbins_y),
  nbins_y = nbins;
end
edges_y = get_edges(fn_y,edges_y,nbins_y,lim_y,lim_prctile_y,outputfun_y,binmode_y);
centers_y = (edges_y(1:end-1)+edges_y(2:end))/2;

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

% flyidx i corresponds to fly flies(i)
fly2idx = sparse(ones(1,nflies),flies,1:nflies,1,obj.nflies);
% expidx i corresponds to experiment ns(i)
n2idx = sparse(ones(1,nexpdirs),ns,1:nexpdirs,1,obj.nexpdirs);

% set expdirs, movie2flies, fly2movies to be relative only to the flies,
% exps selected
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

%% per-fly histograms

% these won't all be used for all averaging conditions, but it should be
% quick and will speed up jackknifing, so compute them
Zperfly = sum(sum(countsperfly,2),3);
fracperfly = bsxfun(@rdivide,countsperfly,Zperfly);

%% compute frac using all the data

% compute final histogram using all the data
frac = compute_frac(1:nflies,1:nexpdirs);

res.frac = frac;

%% jackknife to get an estimate of standard deviation

if docomputestd && ~strcmpi(jackknife,'none'),

  jackknife_frac_mean = zeros([nbins_y,nbins_x]);
  jackknife_frac_std = zeros([nbins_y,nbins_x]);
  nsets = zeros([nbins_y,nbins_x]);
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

% store for returning
res.jackknife_frac_std = jackknife_frac_std;

%% jackknife to get an estimate of standard error
if docomputestderr && ~strcmpi(jackknife,'none'),

  jackknife_frac_mean = zeros([nbins_y,nbins_x]);
  jackknife_frac_stderr = zeros([nbins_y,nbins_x]);
  nsets = zeros([nbins_y,nbins_x]);
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

% store for returning
res.jackknife_frac_stderr = jackknife_frac_stderr;
    
%% plot

res.plotparams = struct;

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
    else
      % handle input, but figure does not exist, so make it
      figure(hfig(figi_main));
      if ~isempty(figpos),
        set(hfig(figi_main),'Position',figpos);
      end
    end
    % create the axes
    hax(axi_main) = createsubplots(1,nax_main,.05,hfig(figi_main));
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

  % TODO: move this
  % user data for main plot
  plotinfo = struct;
  plotinfo.expdirs = expdirs;
  plotinfo.flies = flies;
  plotinfo.averaging = averaging;
  plotinfo.jackknife = jackknife;
  set(hax(axi_hist),'UserData',plotinfo);

  % for clim
  minv = inf;
  maxv = -inf;
  
  %% plot the main histogram
  if fraclogscale,
    im = log(frac);
  else
    im = frac;
  end
  hhist = imagesc(edges_x([1,nbins_x+1]),...
    edges_y([1,nbins_y+1]),im,'parent',hax(axi_hist),...
    'tag','mainhist','userdata',plotinfo);
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
    
    errinfo = struct;
    errinfo.jackknife = jackknife;
    if strcmpi(jackknife,'perfly'),
      errinfo.n = nflies;
    elseif strcmpi(jackknife,'perexp'),
      errinfo.n = nexpdirs;
    end
    
    if fraclogscale,
      im = log(frac-err);
    else
      im = frac-err;
    end
    herr_minus = imagesc(edges_x([1,nbins_x+1]),...
      edges_y([1,nbins_y+1]),im,'parent',hax(axi_minuserr),...
      'tag',sprintf('minuserr_%s',ploterrorbars),'userdata',errinfo);
    % update limits
    minv = min(minv,min(im(~isinf(im))));

    if fraclogscale,
      im = log(frac+err);
    else
      im = frac+err;
    end
    herr_plus = imagesc(edges_x([1,nbins_x+1]),...
      edges_y([1,nbins_y+1]),im,'parent',hax(axi_pluserr),...
      'tag',sprintf('pluserr_%s',ploterrorbars),'userdata',errinfo);

    maxv = max(maxv,max(im(~isinf(im))));

    res.plotparams.herr_minus = herr_minus;
    res.plotparams.herr_plus = herr_plus;
    
    title(hax(axi_minuserr),sprintf('- %s',ploterrorbars));
    title(hax(axi_pluserr),sprintf('+ %s',ploterrorbars));
        
  end

  %% create axes for experiments

  nax_perexp = 0;
  figi_perexp = figi_main;
  axi_off = nax_main;
  if nexpdirs > 1 && plotperexp,
    figi_perexp = figi_main + 1;
    nax_perexp = nexpdirs;
    figpos_perexp = get(hfig(figi_main),'Position');
    figpos_perexp(4) = figpos_perexp(4)*nax_main/nax_perexp;
    if numel(hax) < axi_off + nax_perexp,
      if numel(hfig) < figi_perexp,
        hfig(figi_perexp) = figure;
        set(hfig(figi_perexp),'Position',figpos_perexp);
      elseif ishandle(hfig(figi_perexp)),
        clf(hfig(figi_perexp));
      else
        figure(hfig(figi_perexp));
        set(hfig(figi_perexp),'Position',figpos_perexp);
      end
      hax = [hax,createsubplots(1,nax_perexp,.05,hfig(figi_perexp))];
    else
      hfig(figi_perexp) = get(hax(axi_off+1),'Parent');
    end
    for axi = axis_off+1:axi_off+nax_perexp,
      cla(hax(axi));
      %hold(hax(axi),'on');
    end
  end

  %% plot per-experiment histogram
  
  if plotperexp && nexpdirs > 1,
    hperexphist = zeros(1,nexpdirs);
    for n = 1:nexpdirs,
      expinfo = struct('n',ns(n),'flies',flies(movie2flies{n}),'expdir',expdirs{n},'ndataperfly',Zperfly(movie2flies{n}));
      fracperexpcurr = compute_frac(1:nflies,n);
      if fraclogscale,
        fracperexpcurr = log(fracperexpcurr);
      else
        %fracperexpcurr = fracperexpcurr;
      end

      % update limits
      minv = min(minv,min(fracperexpcurr(~isinf(fracperexpcurr))));
      maxv = max(maxv,max(fracperexpcurr(~isinf(fracperexpcurr))));
      hperexphist(n) = ...
        imagesc(edges_x([1,nbins_x+1]),...
        edges_y([1,nbins_y+1]),fracperexpcurr,'parent',hax(axi_off+n),...
        'tag',sprintf('exp%d',n),'userdata',expinfo);
      title(hax(axi_off+n),expdirs{n},'interpreter','none');
    end
    res.plotparams.hperexphist = hperexphist;
    
  end
  
  %% create axes for flies
  
  if nflies > 1 && plotperfly,
    figi_perfly = figi_perexp + 1;
    axi_off = axi_off + nax_perexp;
    nax_perfly = nflies;
    figpos_perfly = get(hfig(figi_main),'Position');
    figpos_perfly(4) = figpos_perfly(4)*nax_main/nax_perfly;
    if numel(hax) < axi_off + nax_perfly,
      if numel(hfig) < figi_perfly,
        hfig(figi_perfly) = figure;
        set(hfig(figi_perfly),'Position',figpos_perfly);
      elseif ishandle(hfig(figi_perfly)),
        clf(hfig(figi_perfly));
      else
        figure(hfig(figi_perfly));
        set(hfig(figi_perfly),'Position',figpos_perfly);
      end
      hax = [hax,createsubplots(1,nax_perfly,[[.01,.01];[.1,.1]],hfig(figi_perfly))];
    else
      hfig(figi_perfly) = get(hax(axi_off+1),'Parent');
    end
    for axi = axi_off+1:axi_off+nax_perfly,
      cla(hax(axi));
      %hold(hax(axi),'on');
    end
  end
      
  res.plotparams.hax = hax;
  res.plotparams.hfig = hfig;

  %% plot the per-fly histograms
  if plotperfly && nflies > 1,
    hperflyhist = zeros(1,nflies);
    for fly = 1:nflies,
      flyglobal = flies(fly);
      nglobal = obj.fly2movie(flyglobal);
      expdir = obj.expdir_bases{nglobal};
      flyofmovie = find(flyglobal == obj.movie2flies{nglobal},1);
      flyinfo = struct('fly',flyglobal,'n',nglobal,'expdir',expdir,'ndata',Zperfly(fly),'flyofmovie',flyofmovie);
      
      if fraclogscale,
        im = log(fracperfly(fly,:,:));
      else
        im = fracperfly(fly,:,:);
      end

      hperflyhist(fly) = ...
        imagesc(edges_x([1,nbins_x+1]),...
        edges_y([1,nbins_y+1]),permute(im,[2,3,1]),...
        'parent',hax(axi_off+fly),...
        'tag',sprintf('fly%d',fly),'userdata',flyinfo);
      title(hax(axi_off+fly),sprintf('Fly %d of %s',flyofmovie,expdir),'interpreter','none');
      axis(hax(axi_off+fly),'off');
    end
    
    % update limits
    if fraclogscale,
      minv = min(minv,log(min(fracperfly(fracperfly>0))));
      maxv = max(maxv,log(max(fracperfly(fracperfly>0))));
    else
      minv = min(minv,min(fracperfly(:)));
      maxv = max(maxv,max(fracperfly(:)));
    end


    res.plotparams.hperflyhist = hperflyhist;
    
  end
  
  %% set axis limits

  if isnan(clim(1)),
    clim(1) = minv;
  end
  if isnan(clim(2)),
    clim(2) = maxv;
  end
  % make sure the y lims are positive
  if fraclogscale,
    clim = max(clim,log(min(frac(frac>0)))/2);
  end
  for axi = 1:numel(hax),
    set(hax(axi),'clim',clim);
  end
  linkaxes(hax);
  linkprop(hax,'clim');
  
  res.plotparams.clim = clim;

  %% labels

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
  for axi = 1:numel(hax),
    if ~isempty(axesstyleparams),
      set(hax(axi),axesstyleparams{:});
    end
    hxlabel(axi) = xlabel(hax(i),s_x,'interpreter','none');
    hylabel(axi) = ylabel(hax(i),s_y,'interpreter','none');
  end
  plottitle = sprintf('%s vs %s',s_x,s_y);
  set(hfig(figi_main),'name',plottitle);
  if plotperexp && nexpdirs > 1,
    set(hfig(figi_perexp),'name',[plottitle,' per exp']);
  end
  if plotperfly && nflies > 1,
    set(hfig(figi_perfly),'name',[plottitle,' per fly']);
  end
  
  res.hxlabel = hxlabel;
  res.hylabel = hylabel;
    
end

%% compute_frac

function my_frac = compute_frac(my_flies,my_ns)

  % take the intersection of the specified flies and experiments
  my_allflies_perexp = [movie2flies{my_ns}];
  my_flies = intersect(my_flies,my_allflies_perexp);
  my_ns = unique(fly2movie(my_flies));
  
  switch lower(averaging),
    case 'allexps_allflies',
      
      % treat all frames of data the same: just add up all the counts
      my_counts = sum(countsperfly(my_flies,:,:),1);
      my_frac = my_counts / sum(my_counts(:));
      
    case 'allexps_perfly',
      
      % get per-fly fracs, but treat all flies in all experiments the same
      my_frac = nanmean(fracperfly(my_flies,:,:),1);
      
    case 'perexp_allflies',
      
      % get per-exp fracs, but treat all frames within an experiment the
      % same
      
      % loop over experiments
      my_frac = zeros([1,nbins_y,nbins_x]);
      my_nexps = zeros([1,nbins_y,nbins_x]);
      for my_n = my_ns(:)',
        % which flies in this experiment
        my_flies_curr = intersect(my_flies,movie2flies{my_n});
        % frac for these flies
        my_frac_curr = sum(countsperfly(my_flies_curr,:,:),1);
        my_frac_curr = my_frac_curr / sum(my_frac_curr(:));
        % only include good data
        my_isdata = ~isnan(my_frac_curr);
        my_frac(my_isdata) = my_frac(my_isdata) + my_frac_curr(my_isdata);
        my_nexps(my_isdata) = my_nexps(my_isdata) + 1;
      end
      % normalize
      my_frac = my_frac ./ my_nexps;
      
    case 'perexp_perfly',
      
      % loop over experiments
      my_frac = zeros([1,nbins_y,nbins_x]);
      my_nexps = zeros([1,nbins_y,nbins_x]);
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

  my_frac = permute(my_frac,[2,3,1]);

end

function [edges,nbins] = get_edges(fn,edges,nbins,lim,lim_prctile,outputfun,binmode)
  
% set edges if not input
if ~isempty(edges),
  nbins = length(edges);
  return;
end
if any(isnan(lim)),
  minx = inf;
  maxx = -inf;
  % for now, don't worry about memory, store all the data
  useprctile = any(~isnan(lim_prctile));
  if useprctile,
    xall = [];
  end
  for my_fly = flies,
    if ~isempty(conditions),
      my_x = obj.trx(my_fly).(fn)(conditions(obj.trx(my_fly)));
    else
      my_x = obj.trx(my_fly).(fn);
    end
    if ~isempty(outputfun),
      my_x = outputfun(my_x);
    end
    badidx = isinf(my_x) | isnan(my_x);
    minx = min(minx,min(my_x(~badidx)));
    maxx = max(maxx,max(my_x(~badidx)));
    if useprctile,
      xall = [xall,my_x]; %#ok<AGROW>
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

end

end