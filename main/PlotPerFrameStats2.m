function handles = PlotPerFrameStats2(stats_perframefeatures,singlestats,combinedstats,controlstats,basename,varargin)

% parse plotting parameters
[hfig,visible,position,axposition,allfields,...
  allframeconditions,allflyconditions,...
  stattype,weighttype,nframestotal,plotstderr,...
  datanames] = ...
  myparse(varargin,'hfig',1,...
  'visible','off',...
  'position',[1 1 2000 720],...
  'axposition',[.025,.3,.95,.65],...
  'fields',{},...
  'frameconditions',{},...
  'flyconditions',{},...
  'stattype','flymeans',...
  'weighttype','fracframesfly',...
  'nframestotal',[],...
  'plotstderr',true,...
  'datanames',{});

% set up figure
if ~ishandle(hfig),
  figure(hfig);
else
  clf(hfig);
end
set(hfig,'Visible',visible,'Units','pixels','Position',position);
hax = axes('Position',axposition,'Parent',hfig);

% this will be one if we are not comparing multiple results
ndata = numel(combinedstats);

%% fns corresponding to fields

if isempty(allfields),
  allfields = unique({stats_perframefeatures.field});
end
nfields = numel(allfields);

fns = {};
flyconditions = {};
frameconditions = {};
fields = {};
% only allow conditions in allframeconditions, allflyconditions
goodidx = true(1,numel(stats_perframefeatures));
if ~isempty(allframeconditions),
  goodidx = ismember({stats_perframefeatures.framecondition},allframeconditions);
end
if ~isempty(allflyconditions),
  goodidx = goodidx & ismember({stats_perframefeatures.flycondition},allflyconditions);
end

for fieldi = 1:nfields,
  field = allfields{fieldi};
  idx = find(strcmp(field,{stats_perframefeatures.field}) & goodidx);
  for i = idx,
    fns{end+1} = sprintf('%s_fly%s_frame%s',field,stats_perframefeatures(i).flycondition,...
      stats_perframefeatures(i).framecondition); %#ok<AGROW>
    % truncate name
    if numel(fns{end}) > 63,
      fns{end} = fns{end}(1:63);
    end
    fields{end+1} = field; %#ok<AGROW>
  end
  flyconditions = [flyconditions,{stats_perframefeatures(idx).flycondition}]; %#ok<AGROW>
  frameconditions = [frameconditions,{stats_perframefeatures(idx).framecondition}]; %#ok<AGROW>
end

nfns = numel(fns);

%% z-score combinedstats

if iscell(combinedstats),
  nprctiles = numel(combinedstats{1}.(fns{1}).meanprctiles);
else
  nprctiles = numel(combinedstats(1).(fns{1}).meanprctiles);
end

meanmeans = nan(ndata,nfns);
stdmeans = nan(ndata,nfns);
stderrmean = nan(ndata,nfns);
meanprctiles = nan(nprctiles,nfns,ndata);
stdprctiles = nan(nprctiles,nfns,ndata);
stderrprctiles = nan(nprctiles,nfns,ndata);
nanalyzed = nan(ndata,nfns);
stderrmean_control = nan(1,nfns);

for fni = 1:nfns,

  fn = fns{fni};
  
  % control data mean and standard deviation
  if isempty(controlstats),
    mu = 0;
    sig = 1;
    stderrmean_control(fni) = 1;
  else
    if strcmpi(stattype,'flymeans'),
      datacurrmean = controlstats.meanstatsperfly.(fn);
      datacurr = controlstats.statsperfly.(fn);
    else
      datacurrmean = controlstats.meanstatsperexp.(fn);
      datacurr = controlstats.statsperexp.(fn);
    end
    mu = datacurrmean.meanmean;
    sig = datacurrmean.stdmean;
    goodidx = ~isnan(datacurr.Z) & datacurr.Z > 0;
    if strcmpi(stattype,'flymeans'),
      switch lower(weighttype),
        case 'nframesfly'
          if isempty(nframestotal),
            error('nframestotal required for weighttype %s',weighttype);
          end
          ncontrol = sum(datacurr.Zfly(goodidx)) / nframestotal;
        case 'fracframesfly',
          ncontrol = sum(datacurr.fracframesanalyzed(goodidx));
        case 'nframesanalyzed',
          ncontrol = nan;
        case 'uniform',
          ncontrol = nnz(goodidx);
      end
    else
      ncontrol = nnz(goodidx);
    end
    stderrmean_control(fni) = 1/sqrt(ncontrol);
  end

  for datai = 1:ndata,
    if iscell(combinedstats),
      datacurr = combinedstats{datai}.(fn);
    else
      datacurr = combinedstats(datai).(fn);
    end
      
    meanmeans(datai,fni) = (datacurr.meanmean - mu) / sig;
    stdmeans(datai,fni) = datacurr.stdmean / sig;
    meanprctiles(:,fni,datai) = (datacurr.meanprctiles - mu) / sig;
    stdprctiles(:,fni,datai) = datacurr.stdprctiles / sig;
    
    if iscell(combinedstats),
      datacurr = singlestats{datai}.(fn);
    else
      datacurr = singlestats(datai).(fn);
    end
    goodidx = ~isnan(datacurr.Z) & datacurr.Z > 0;
    
    if strcmpi(stattype,'flymeans'),
      switch lower(weighttype),
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
    
    nanalyzed(datai,fni) = n;
    stderrmean(datai,fni) = stdmeans(datai,fni) / sqrt(n);
    stderrprctiles(:,fni,datai) = stdprctiles(:,fni,datai) / sqrt(n);
  
  end
end

%% colors

if ndata == 1,
  colors = jet(nfns)*.7;
  [~,order] = sort(min(abs(meanmeans+stderrmean),abs(meanmeans-stderrmean)));
  [~,order] = sort(order);
  colors = colors(order,:);
elseif ndata > 7,
  colors = jet(ndata)*.7;
else
  colors = lines(ndata);
end

%% plot std/stderr of control data

hold(hax,'on');
if ~isempty(controlstats) && ndata == 1,
  if plotstderr,
    d = stderrmean_control;
  else
    d = 1;
  end
  hcontrol(1) = plot(hax,1:nfns,d,':','color',[.5,.5,.5]);
  hcontrol(2) = plot(hax,1:nfns,-d,':','color',[.5,.5,.5]);
else
  hcontrol = [];
end

%% plot std/stderr, mean of means

if plotstderr,
  d = stderrmean;
else
  d = stdmeans;
end

if ndata == 1,
  h = nan(1,nfns);
  for i = 1:nfns,
    h(i) = errorbar(hax,i,meanmeans(datai,i),d(datai,i),'o','color',colors(i,:),'markerfacecolor',colors(i,:));
  end
else
  h = nan(1,ndata);
  for datai = 1:ndata,
    plot(repmat(1:nfns,[2,1]),bsxfun(@plus,meanmeans(datai,:),[d(datai,:);-d(datai,:)]),'-','color',colors(datai,:));
    h(datai) = plot(1:nfns,meanmeans(datai,:),'o','color',colors(datai,:),'markerfacecolor',colors(datai,:));
  end
end

%% set axis limits

set(hax,'XLim',[0,nfns+1]);

maxy = max(meanmeans(:)+d(:));
miny = min(meanmeans(:)-d(:));
ylim = [miny-(maxy-miny)*.01,maxy+(maxy-miny)*.01];
set(hax,'YLim',ylim);

%% tick labels

allflyconditions = unique(flyconditions);
allframeconditions = unique(frameconditions);
if numel(allflyconditions) == 1 && numel(allframeconditions) == 1,
  typestrs = fields;
elseif numel(allflyconditions) == 1,
  typestrs = cell(1,nfns);
  for i = 1:nfns,
    typestrs{i} = sprintf('%s,%s',fields{i},frameconditions{i});
  end
elseif numel(allframeconditions) == 1,
  typestrs = cell(1,nfns);
  for i = 1:nfns,
    typestrs{i} = sprintf('%s,%s',fields{i},flyconditions{i});
  end
else
  typestrs = cell(1,nfns);
  for i = 1:nfns,
    typestrs{i} = sprintf('%s, %s, %s',fields{i},flyconditions{i},frameconditions{i});
  end
end
set(hax,'XTick',1:nfns,'XTickLabel',typestrs);

%% legend

%xlabel(hax,'Experiment','Interpreter','none');
if isempty(controlstats),
  hy = ylabel(hax,'Behavior statistic','Interpreter','None');
else
  hy = ylabel(hax,'Control stds','Interpreter','None');
end
if plotstderr,
  if strcmpi(stattype,'flymeans'),
    hti = title(hax,sprintf('Mean, stderr of per-fly means for %s',basename),'Interpreter','none');
  else
    hti = title(hax,sprintf('Mean, stderr of per-experiment means for %s',basename),'Interpreter','none');
  end
else
  if strcmpi(stattype,'flymeans'),
    hti = title(hax,sprintf('Mean, std of per-fly means for %s',basename),'Interpreter','none');
  else
    hti = title(hax,sprintf('Mean, std of per-experiment means for %s',basename),'Interpreter','none');
  end
end

if ndata > 1 && ~isempty(datanames),
  hle = legend(hax,h,datanames);
  set(hle,'Interpreter','none');
else
  hle = [];
end

%% rotate ticks
hx = rotateticklabel(hax,90);
% make sure the ticks don't overlap the x-axis
ex = get(hx(1),'Extent');
y1 = ex(2)+ex(4);
offy = max(0,y1-miny);
for i = 1:numel(hx),
  pos = get(hx(i),'Position');
  pos(2) = pos(2) - offy;
  set(hx(i),'Position',pos);
  if ndata == 1,
    set(hx(i),'color',colors(i,:));
  end
end
set(hx,'Interpreter','none');
set(hax,'box','off');
if ndata > 1,
  set(hax,'xgrid','on');
end
%% outputs

handles = struct;
handles.hfig = hfig;
handles.hax = hax;
handles.hcontrol = hcontrol;
handles.h = h;
handles.hy = hy;
handles.hti = hti;
handles.hx = hx;
handles.hle = hle;