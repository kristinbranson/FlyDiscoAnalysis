function hfig = ComparePerFrameStats(expdirs,field,varargin)

%% parse arguments
[settingsdir,analysis_protocol,datalocparamsfilestr,...
  hfig,position,axposition,sortby,offx,colors] = ...
  myparse(varargin,...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'analysis_protocol','current',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'hfig',1,...
  'position',[1 1 1500 500],...
  'axposition',[.1,.1,.85,.8],...
  'sortby','none',...
  'offx',.2,...
  'colors',zeros(0,3));

nexpdirs = numel(expdirs);
basenames = cell(1,nexpdirs);
for i = 1:nexpdirs,
  [~,basenames{i}] = fileparts(expdirs{i});
end

%% file locations
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

%% read plotting parameters
% statsplotparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statsplotparamsfilestr);
% stats_plot_params = ReadParams(statsplotparamsfile);
statsperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statsperframefeaturesfilestr);
stats_perframefeatures = ReadStatsPerFrameFeatures(statsperframefeaturesfile);
histplotparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histplotparamsfilestr);
hist_plot_params = ReadParams(histplotparamsfile);

%% fns corresponding to field

idx = find(strcmp(field,{stats_perframefeatures.field}));
nfns = numel(idx);
fns = cell(1,nfns);
for ii = 1:nfns,
  i = idx(ii);
  fns{ii} = sprintf('%s_fly%s_frame%s',field,stats_perframefeatures(i).flycondition,...
    stats_perframefeatures(i).framecondition);
end
flyconditions = {stats_perframefeatures(idx).flycondition};
frameconditions = {stats_perframefeatures(idx).framecondition};

%% set up figure
black = ones(1,3)/255;

if ~ishandle(hfig),
  figure(hfig);
else
  clf(hfig);
end
set(hfig,'Position',position);
hax = axes('Position',axposition,'Parent',hfig,'XColor',black,'YColor',black);
hold(hax,'on');

%% choose colors for experiments


colors = cat(1,colors,jet(nexpdirs-size(colors,1))*.7);

%% load data for all experiments

% get nprctiles
expdiri = 1;
expdir = expdirs{expdiri};
statsmatfile = fullfile(expdir,dataloc_params.statsperframematfilestr);
load(statsmatfile,'statsperexp');
nprctiles = numel(statsperexp.(fns{1}).meanprctiles);

meanmeans = nan(nfns,nexpdirs);
stdmeans = nan(nfns,nexpdirs);
stderrmean = nan(nfns,nexpdirs);
meanprctiles = nan(nprctiles,nexpdirs,nfns);
stdprctiles = nan(nprctiles,nexpdirs,nfns);
stderrprctiles = nan(nprctiles,nexpdirs,nfns);
nfliesanalyzed = nan(nexpdirs,nfns);

for expdiri = 1:nexpdirs,
  
  expdir = expdirs{expdiri};
  
  % load data
  statsmatfile = fullfile(expdir,dataloc_params.statsperframematfilestr);
  load(statsmatfile,'statsperexp','statsperfly');

  % store
  for fni = 1:nfns,
    fn = fns{fni};
    meanmeans(fni,expdiri) = statsperexp.(fn).meanmean;
    stdmeans(fni,expdiri) = statsperexp.(fn).stdmean;
    meanprctiles(:,expdiri,fni) = statsperexp.(fn).meanprctiles;
    stdprctiles(:,expdiri,fni) = statsperexp.(fn).stdprctiles;
    
    goodidx = ~isnan(statsperfly.(fn).Z);
    nflies = sum(statsperfly.(fn).fracframesanalyzed(goodidx));
    nfliesanalyzed(expdiri,fni) = nflies;
    stderrmean(fni,expdiri) = stdmeans(fni,expdiri) / sqrt(nflies);
    stderrprctiles(:,expdiri,fni) = stdprctiles(:,expdiri,fni) / sqrt(nflies);
    
  end
  
end

%% sort experiments

switch lower(sortby),
  
  case 'alphabetical',
    [~,order] = sort(basenames);
  case 'largest',
    [~,fnisort] = max(nansum(nfliesanalyzed,1));
    [~,order] = sort(meanmeans(fnisort,:));
  otherwise
    order = 1:nexpdirs;
end

expdirs = expdirs(order);
basenames = basenames(order);
meanmeans = meanmeans(:,order);
stdmeans = stdmeans(:,order);
stderrmean = stderrmean(:,order);
meanprctiles = meanprctiles(:,order,:);
stdprctiles = stdprctiles(:,order,:);
stderrprctiles = stderrprctiles(:,order,:);
nfliesanalyzed = nfliesanalyzed(order,:);

%% plot

% plot stderr of mean means
hstd = nan(1,nexpdirs);

for expdiri = 1:nexpdirs,
  if nexpdirs == 1,
    offcurr = 0;
  else
    offcurr = (2*(expdiri-(nexpdirs+1)/2)/(nexpdirs-1))*offx;
  end
  hstd(expdiri) = errorbar(hax,(1:nfns)+offcurr,meanmeans(:,expdiri),stderrmean(:,expdiri),'.','color',colors(expdiri,:));
end

% plot the mean means
hmean = nan(1,nexpdirs);
for expdiri = 1:nexpdirs,
  if nexpdirs == 1,
    offcurr = 0;
  else
    offcurr = (2*(expdiri-(nexpdirs+1)/2)/(nexpdirs-1))*offx;
  end
  hmean(expdiri) = plot(hax,(1:nfns)+offcurr,meanmeans(:,expdiri),'o:','color',colors(expdiri,:),'markerfacecolor',colors(expdiri,:));
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
if numel(allflyconditions) == 1,
  typestrs = frameconditions;
elseif numel(allframeconditions) == 1,
  typestrs = flyconditions;
else
  typestrs = cell(1,nfns);
  for i = 1:nfns,
    typestrs{i} = sprintf('%s, %s',flyconditions{i},frameconditions{i});
  end
end
set(hax,'XTick',1:nfns,'XTickLabel',typestrs);

%% legend

legend(hmean,basenames,'Location','EastOutside','Parent',hfig,'Interpreter','none');
%xlabel(hax,'Experiment','Interpreter','none');
ylabel(hax,field,'Interpreter','None');
title(hax,sprintf('Mean, stderr of mean for %s',field),'Interpreter','none');

