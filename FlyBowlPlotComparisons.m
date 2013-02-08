function FlyBowlPlotComparisons(rootdatadir,datadirs,comparison_name,varargin)

[analysis_protocol,settingsdir,datalocparamsfilestr,controlstats,visible,...
  meanweighttype,stdweighttype,nframestotal,fly_plotstderr,exp_plotstderr] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'controlstats',[],...
  'visible','off',...
  'meanweighttype','nframesfly',...
  'stdweighttype','fracframesfly',...
  'nframestotal',[],...
  'fly_plotstderr',true,...
  'exp_plotstderr',true);

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
statsperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statsperframefeaturesfilestr);
stats_perframefeatures = ReadStatsPerFrameFeatures2(statsperframefeaturesfile);
statframeconditionsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statframeconditionfilestr);
frameconditiondict = ReadParams(statframeconditionsfile);
%statflyconditionsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statflyconditionfilestr);
%flyconditiondict = ReadParams(statflyconditionsfile);
%statsparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statsparamsfilestr);
%stats_params = ReadParams(statsparamsfile);

histplotparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histplotparamsfilestr);
hist_plot_params = ReadParams(histplotparamsfile);
histperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframefeaturesfilestr);
hist_perframefeatures = ReadHistPerFrameFeatures2(histperframefeaturesfile);
histperframebinsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframebinsfilestr);
load(histperframebinsfile,'bins');

%% load in data

ncomparisons = numel(datadirs);

allmeanstatsperfly = cell(1,ncomparisons);
allstatsperfly = cell(1,ncomparisons);
allmeanstatsperexp = cell(1,ncomparisons);
allstatsperexp = cell(1,ncomparisons);
for i = 1:ncomparisons,
  datacurr = load(fullfile(rootdatadir,datadirs{i},dataloc_params.statsperframematfilestr));
  allstatsperfly{i} = datacurr.statsperfly;
  allmeanstatsperfly{i} = datacurr.meanstatsperfly;
  allstatsperexp{i} = datacurr.statsperexp;
  allmeanstatsperexp{i} = datacurr.meanstatsperexp;
end

%% plot all stats together

fprintf('Plotting overview stats...\n');
handles = PlotPerFrameStats2(stats_perframefeatures,allstatsperfly,allmeanstatsperfly,controlstats,comparison_name,...
  'visible',visible,...
  'stattype','flymeans',...
  'weighttype','fracframesfly',...
  'plotstderr',fly_plotstderr,...
  'datanames',datadirs);

%handles = PlotPerFrameStatsComparison(stats_perframefeatures,allstatsperfly,allmeanstatsperfly,controlstats,datadirs);
filename = 'stats';
exts = {'pdf','jpeg'};
if strcmpi(comparison_name,'all'),
  outdir = fullfile(rootdatadir,'analysis_plots_new');
else
  outdir = fullfile(rootdatadir,comparison_name);
end
if ~exist(outdir,'dir'),
  mkdir(outdir);
end
for i = 1:numel(exts),
  savefig([filename,'.',exts{i}],handles.hfig,exts{i});
  unix(sprintf('mv %s %s',[filename,'.',exts{i}],outdir));
end
set(handles.hfig,'Units','pixels');

handles = PlotPerFrameStats2(stats_perframefeatures,allstatsperexp,allmeanstatsperexp,controlstats,comparison_name,...
  'visible',visible,...
  'stattype','expmeans',...
  'plotstderr',exp_plotstderr,...
  'datanames',datadirs);

%handles = PlotPerFrameStatsComparisonExp(stats_perframefeatures,allstatsperexp,allmeanstatsperexp,controlstats,datadirs,'plotstderr',exp_plotstderr);
filename = 'stats_exp';
for i = 1:numel(exts),
  savefig([filename,'.',exts{i}],handles.hfig,exts{i});
  unix(sprintf('mv %s %s',[filename,'.',exts{i}],outdir));
end

%% plot all histograms together

allmeanhistperfly = cell(1,ncomparisons);
allhistperfly = cell(1,ncomparisons);
allmeanhistperexp = cell(1,ncomparisons);
allhistperexp = cell(1,ncomparisons);
for i = 1:ncomparisons,
  datacurr = load(fullfile(rootdatadir,datadirs{i},dataloc_params.histperframematfilestr));
  allhistperfly{i} = datacurr.histperfly;
  allmeanhistperfly{i} = datacurr.meanhistperfly;
  allhistperexp{i} = datacurr.histperexp;
  allmeanhistperexp{i} = datacurr.meanhistperexp;
end

%exts = {'jpeg','pdf'};
%outdir = fullfile(rootdatadir,'analysis_plots');

for typei = 1:numel(hist_perframefeatures),
  
  fprintf('Plotting histograms for field=%s, fly=%s, frame=%s...\n',...
    hist_perframefeatures(typei).field,...
    hist_perframefeatures(typei).flycondition,...
    hist_perframefeatures(typei).framecondition);
  
  if strcmp(hist_perframefeatures(typei).field,'duration'),
    binfn = ['duration_',hist_perframefeatures(typei).framecondition];    
%     frameconditionparams = DecodeConditions(hist_perframefeatures(typei).framecondition,frameconditiondict);
%     m = regexp(frameconditionparams(1:2:end),'^[^_]+_(.+)_labels$','once','tokens');
%     tmp = find(~cellfun(@isempty,m),1);
%     binfn = ['duration','_',m{tmp}{1}];
  else
    binfn = hist_perframefeatures(typei).field;
  end
  
%   handles = PlotPerFrameHistsComparison(binfn,hist_perframefeatures(typei).field,...
%     hist_perframefeatures(typei).flycondition,...
%     hist_perframefeatures(typei).framecondition,...
%     allmeanhistperfly,allhistperfly,...
%     bins.(binfn),...
%     hist_plot_params,...
%     datadirs,...
%     'visible',visible);
%     
  handles = PlotPerFrameHists3(binfn,hist_perframefeatures(typei).field,...
    typei,hist_perframefeatures,...
    allmeanhistperfly,allhistperfly,...
    bins.(binfn),hist_plot_params,...
    comparison_name,...
    'datanames',datadirs,'stattype','flymeans',...
    'meanweighttype',meanweighttype,'stdweighttype',stdweighttype,...
    'nframestotal',nframestotal,...
    'plotstderr',fly_plotstderr,...
    'visible',visible);
  
  drawnow;
  filename = sprintf('hist_%s_fly%s_frame%s',...
    hist_perframefeatures(typei).field,...
    hist_perframefeatures(typei).flycondition,...
    hist_perframefeatures(typei).framecondition);
  for i = 1:numel(exts),
    savefig([filename,'.',exts{i}],handles.hfig,exts{i});
    unix(sprintf('mv %s %s',[filename,'.',exts{i}],outdir));
  end
  set(handles.hfig,'Units','pixels');

%   handles = PlotPerFrameHistsComparisonExp(binfn,hist_perframefeatures(typei).field,...
%     hist_perframefeatures(typei).flycondition,...
%     hist_perframefeatures(typei).framecondition,...
%     allmeanhistperexp,allhistperexp,...
%     bins.(binfn),...
%     hist_plot_params,...
%     datadirs,...
%     'visible',visible,...
%     'plotstderr',exp_plotstderr);
  
  handles = PlotPerFrameHists3(binfn,hist_perframefeatures(typei).field,...
    typei,hist_perframefeatures,...
    allmeanhistperexp,allhistperexp,...
    bins.(binfn),hist_plot_params,...
    comparison_name,...
    'datanames',datadirs,'stattype','expmeans',...
    'plotstderr',exp_plotstderr,...
    'visible',visible);

  drawnow;
  
  filename = sprintf('hist_%s_fly%s_frame%s_exp',...
    hist_perframefeatures(typei).field,...
    hist_perframefeatures(typei).flycondition,...
    hist_perframefeatures(typei).framecondition);
  for i = 1:numel(exts),
    savefig([filename,'.',exts{i}],handles.hfig,exts{i});
    unix(sprintf('mv %s %s',[filename,'.',exts{i}],outdir));
  end
  set(handles.hfig,'Units','pixels');


end
% 
% handles = PlotPerFrameHists(field,hist_perframefeatures,...
%     meanhistperfly,histperfly,...
%     bins.(field),hist_plot_params,plottitle,...
%     'visible',visible,...
%     'hax',hax);

%% end loop over comparisons