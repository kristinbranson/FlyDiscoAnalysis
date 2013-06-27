function handles = PlotPerFrameStats(stats_perframefeatures,statsperfly,statsperexp,controlstats,basename,varargin)

% parse plotting parameters
[hfig,visible,position,axposition,allfields,...
  allframeconditions,allflyconditions] = ...
  myparse(varargin,'hfig',1,...
  'visible','off',...
  'position',[1 1 2000 720],...
  'axposition',[.025,.3,.95,.65],...
  'fields',{},...
  'frameconditions',{},...
  'flyconditions',{});

% set up figure
if ~ishandle(hfig),
  figure(hfig);
end

clf(hfig);
set(hfig,'Units','pixels','Visible',visible,'Position',position);
hax = axes('Position',axposition,'Parent',hfig);

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
    fields{end+1} = field;
  end
  flyconditions = [flyconditions,{stats_perframefeatures(idx).flycondition}]; %#ok<AGROW>
  frameconditions = [frameconditions,{stats_perframefeatures(idx).framecondition}]; %#ok<AGROW>
end

nfns = numel(fns);

%% z-score statsperexp

nprctiles = numel(statsperexp.(fns{1}).meanprctiles);

meanmeans = nan(1,nfns);
stdmeans = nan(1,nfns);
stderrmean = nan(1,nfns);
meanprctiles = nan(nprctiles,nfns);
stdprctiles = nan(nprctiles,nfns);
stderrprctiles = nan(nprctiles,nfns);
nfliesanalyzed = nan(1,nfns);
stderrmean_control = nan(1,nfns);

for fni = 1:nfns,

  fn = fns{fni};
  
  % control data mean and standard deviation
  if isempty(controlstats),
    mu = 0;
    sig = 1;
    stderrmean_control(fni) = 1;
  else
    mu = controlstats.meanstatsperfly.(fn).meanmean;
    sig = controlstats.meanstatsperfly.(fn).stdmean;
    goodidx = ~isnan(controlstats.statsperfly.(fn).Z);
    nfliescontrol = sum(controlstats.statsperfly.(fn).fracframesanalyzed(goodidx));
    stderrmean_control(fni) = 1/sqrt(nfliescontrol);
  end

  meanmeans(fni) = (statsperexp.(fn).meanmean - mu) / sig;
  stdmeans(fni) = statsperexp.(fn).stdmean / sig;
  meanprctiles(:,fni) = (statsperexp.(fn).meanprctiles - mu) / sig;
  stdprctiles(:,fni) = statsperexp.(fn).stdprctiles / sig;
  
  goodidx = ~isnan(statsperfly.(fn).Z);
  nflies = sum(statsperfly.(fn).fracframesanalyzed(goodidx));
  nfliesanalyzed(fni) = nflies;
  stderrmean(fni) = stdmeans(fni) / sqrt(nflies);
  stderrprctiles(:,fni) = stdprctiles(:,fni) / sqrt(nflies);
  
end

%% colors

colors = jet(nfns)*.7;
[tmp,order] = sort(min(abs(meanmeans+stderrmean),abs(meanmeans-stderrmean)));
[tmp,order] = sort(order);
colors = colors(order,:);

%% plot stderr of control data

hold(hax,'on');
if ~isempty(controlstats),
  hcontrol(1) = plot(hax,1:nfns,stderrmean_control,':','color',[.5,.5,.5]);
  hcontrol(2) = plot(hax,1:nfns,-stderrmean_control,':','color',[.5,.5,.5]);
else
  hcontrol = [];
end

%% plot stderr, mean of means

h = nan(1,nfns);
for i = 1:nfns,
  h(i) = plot(hax,i,meanmeans(i),'o','color',colors(i,:),'markerfacecolor',colors(i,:));
  h1(i) = plot(hax,[i,i],meanmeans(i)+stderrmean(i)*[-1,1],'-','color',colors(i,:));
end

%% set axis limits

set(hax,'XLim',[0,nfns+1]);

maxy = max(meanmeans(:)+stderrmean(:));
miny = min(meanmeans(:)-stderrmean(:));
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
hti = title(hax,sprintf('Mean, stderr of mean for %s',basename),'Interpreter','none');

%% rotate ticks
hx = rotateticklabel(hax,90);
% make sure the ticks don't overlap the x-axis
ex = get(hx(1),'Extent');
y1 = ex(2)+ex(4);
offy = max(0,y1-miny);
for i = 1:numel(hx),
  pos = get(hx(i),'Position');
  pos(2) = pos(2) - offy;
  set(hx(i),'Position',pos,'color',colors(i,:));
end
set(hx,'Interpreter','none');
set(hax,'box','off');
%% outputs

handles = struct;
handles.hfig = hfig;
handles.hax = hax;
handles.hcontrol = hcontrol;
handles.h = h;
handles.hy = hy;
handles.hti = hti;
handles.hx = hx;
handles.position = position;