% combine data from multiple experiments
function FlyBowlCombineExperiments2(rootdir,outdir,varargin)

special_cases = {'fractime','duration','boutfreq'};

% parse arguments
[settingsdir,analysis_protocol,datalocparamsfilestr,...
  ~,~,visible,controldatadirstr,plottitle,controlstats,controlhist,leftovers] = ...
  myparse_nocheck(varargin,...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'analysis_protocol','current',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'requiredfiles',{'histperframematfilestr','statsperframematfilestr'},...
  'subreadfiles',{},...
  'visible','on',...
  'controldatadirstr','current',...
  'plottitle','',...
  'controlstats',[],...
  'controlhist',[]);


%% get all experiments that satisfy input conditions; currently parsing
% directory structure
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
statsperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statsperframefeaturesfilestr);
stats_perframefeatures = ReadStatsPerFrameFeatures2(statsperframefeaturesfile);
histperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframefeaturesfilestr);
hist_perframefeatures = ReadHistPerFrameFeatures2(histperframefeaturesfile);

%for i = 1:numel(requiredfiles),
%  subreadfiles{end+1} = dataloc_params.(requiredfiles{i}); %#ok<AGROW>
%end
%subreadfiles = unique(subreadfiles);

exp_params = [{'rootdir',rootdir},leftovers];

experiments = SAGEListBowlExperiments(exp_params{:});
nexpdirs = numel(experiments);

if nexpdirs == 0,
  warning('No experiments selected');
  return;
end
expdirs = {experiments.file_system_path};

%% remove bad experiments

isdata = true(1,nexpdirs);
for i = 1:nexpdirs,
  expdir = expdirs{i};
  statsmatfile = fullfile(expdir,dataloc_params.statsperframematfilestr);
  if ~exist(statsmatfile,'file'),
    warning('Stats file %s does not exist',statsmatfile);
    isdata(i) = false;
  end
  histmatfile = fullfile(expdir,dataloc_params.histperframematfilestr);
  if ~exist(histmatfile,'file'),
    warning('Hist file %s does not exist',histmatfile);
    isdata(i) = false;
  end
end
experiments(~isdata) = [];
expdirs(~isdata) = [];
nexpdirs = numel(experiments);

%% load per-fly data

statsperfly = [];
histperfly = [];
statsperexp = [];
histperexp = [];

for i = 1:nexpdirs,
  
  % load the data for the current experiment
  %i=i+1
  fprintf('Loading data for experiment %s\n',expdir);
  expdir = expdirs{i};
  statsmatfile = fullfile(expdir,dataloc_params.statsperframematfilestr);
  statsdata = load(statsmatfile);
  histmatfile = fullfile(expdir,dataloc_params.histperframematfilestr);
  histdata = load(histmatfile);
  
  % merge with data read in so far
  statsperfly = MergeStatsPerFly(statsperfly,statsdata.statsperfly,expdir);
  statsperexp = MergeStatsPerExp(statsperexp,statsdata.statsperexp,expdir);
  %numel(statsperfly.a_mm_flyany_frameany.Z)
  histperfly = MergeHistPerFly(histperfly,histdata.histperfly,expdir);
  histperexp = MergeHistPerExp(histperexp,histdata.histperexp,expdir);

end

%% Combine per-exp data

meanstatsperexp = struct;
meanhistperexp = struct;
fns = fieldnames(statsperexp);
for i = 1:numel(fns),
  fn = fns{i};  
  meanstatsperexp.(fn) = CombinePerExpStats(statsperexp.(fn));
end

fns = fieldnames(histperexp);
for i = 1:numel(fns),
  fn = fns{i};

  % compute per-experiment hist
  meanhistperexp.(fn) = CombinePerExpHists(histperexp.(fn));
  
end

%% Combine per-fly data as in one experiment

meanstatsperfly = struct;
meanhistperfly = struct;

for i = 1:numel(stats_perframefeatures),
  
  field = stats_perframefeatures(i).field;

  % which frames to analyze
  frameconditionname = stats_perframefeatures(i).framecondition;    
  % which flies to analyze
  flyconditionname = stats_perframefeatures(i).flycondition;  
  % minimum number of frames to use this trajectory as a fly
  minZboth = stats_perframefeatures(i).minZboth;
  minZfly = stats_perframefeatures(i).minZfly;
  % unique identifier for this computation, composed of field + conditions
  fn = sprintf('%s_fly%s_frame%s',field,flyconditionname,frameconditionname);
  if numel(fn) > 63,
    fn = fn(1:63);
  end

  % compute per-experiment hist
  meanstatsperfly.(fn) = CombinePerFrameStats2(statsperfly.(fn),minZboth,minZfly);
  
end

for i = 1:numel(hist_perframefeatures),

  field = hist_perframefeatures(i).field;

  % which frames to analyze
  frameconditionname = hist_perframefeatures(i).framecondition;    
  % which flies to analyze
  flyconditionname = hist_perframefeatures(i).flycondition;  
  % minimum number of frames to use this trajectory as a fly
  minZboth = hist_perframefeatures(i).minZboth;
  minZfly = hist_perframefeatures(i).minZfly;
  % unique identifier for this computation, composed of field + conditions
  fn = sprintf('%s_fly%s_frame%s',field,flyconditionname,frameconditionname);
  if numel(fn) > 63,
    fn = fn(1:63);
  end
  
  % compute per-experiment hist
  meanhistperfly.(fn) = CombinePerFrameHists2(histperfly.(fn),minZboth,minZfly);
  
end

%% create output directory
if ~exist(outdir,'file'),
  [success,msg] = mkdir(outdir);
  if ~success,
    error('Could not create directory %s for output:\n%s',outdir,msg);
  end
end

fprintf('Combining %d experiments into directory %s...\n',nexpdirs,outdir);

%% save 

statstxtsavename = fullfile(outdir,dataloc_params.statsperframetxtfilestr);
statsmatsavename = fullfile(outdir,dataloc_params.statsperframematfilestr);

statsdata.statsperfly = statsperfly;
statsdata.statsperexp = statsperexp;
statsdata.meanstatsperfly = meanstatsperfly;
statsdata.meanstatsperexp = meanstatsperexp;
statsdata.exp_params = exp_params;
statsdata.expdirs = expdirs;
statsdata.experiments = experiments;

save(statsmatsavename,'-struct','statsdata');
% save to text file
fns = fieldnames(statsperfly);
for i = 1:numel(fns),
  fn = fns{i};
  SavePerFrameStatsTxtFile(statstxtsavename,fn,statsperfly.(fn),meanstatsperfly.(fn));
end

histtxtsavename = fullfile(outdir,dataloc_params.histperframetxtfilestr);
histmatsavename = fullfile(outdir,dataloc_params.histperframematfilestr);

histdata.histperfly = histperfly;
histdata.histperexp = histperexp;
histdata.meanhistperfly = meanhistperfly;
histdata.meanhistperexp = meanhistperfly;
histdata.exp_params = exp_params;
histdata.expdirs = expdirs;
histdata.experiments = experiments;

save(histmatsavename,'-v7.3','-struct','histdata');
% save to text file
fns = fieldnames(histperfly);
for i = 1:numel(fns),
  fn = fns{i};
  SavePerFrameHistTxtFile(histtxtsavename,fn,histperfly.(fn),meanhistperfly.(fn));
end

%% write experiment info to a file
expinfofile = fullfile(outdir,dataloc_params.combine_experimentsfilestr);
fid = fopen(expinfofile,'w');
if fid < 0,
  error('Could not open file %s for writing',expinfofile);
end

for i = 1:2:numel(exp_params)-1,
  
  s = exp_params{i};
  val = exp_params{i+1};
  if ~iscell(val) && strcmp(s,'not_started'),
    val = {num2str(val)};
  elseif ~iscell(val),
    val = {val};
  end
  
  fprintf(fid,s);
  fprintf(fid,',%s',val{:});
  fprintf(fid,'\n');
end

fprintf(fid,'expdirs');
fprintf(fid,',%s',expdirs{:});
fprintf(fid,'\n');

fclose(fid);

%% create the plot directory if it does not exist
figdir = fullfile(outdir,dataloc_params.figdir);
if ~exist(figdir,'file'),
  [status,msg,~] = mkdir(figdir);
  if ~status,
    error('Could not create the figure directory %s:\n%s',figdir,msg);
  end
end

%% read plotting parameters

histplotparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histplotparamsfilestr);
hist_plot_params = ReadParams(histplotparamsfile);
if isempty(plottitle),
  [~,plottitle] = fileparts(outdir);
end
histperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframefeaturesfilestr);
hist_perframefeatures = ReadHistPerFrameFeatures2(histperframefeaturesfile);
histperframebinsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframebinsfilestr);
load(histperframebinsfile,'bins');
statsperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statsperframefeaturesfilestr);
stats_perframefeatures = ReadStatsPerFrameFeatures2(statsperframefeaturesfile);
statframeconditionsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statframeconditionfilestr);
frameconditiondict = ReadParams(statframeconditionsfile);


%% load control data

if isempty(controlstats) && ~isempty(controldatadirstr),
  controldatadir = fullfile(dataloc_params.pBDPGAL4Ustatsdir,controldatadirstr);
  controlstatsname = fullfile(controldatadir,dataloc_params.statsperframematfilestr);
  controlstats = load(controlstatsname);
end

if isempty(controlhist) && ~isempty(controldatadirstr),
  controldatadir = fullfile(dataloc_params.pBDPGAL4Ustatsdir,controldatadirstr);
  controlhistname = fullfile(controldatadir,dataloc_params.histperframematfilestr);
  controlhist = load(controlhistname);
end

%% plot means, stds

if isstruct(controlstats),
  %stathandles = PlotPerFrameStats(stats_perframefeatures,statsperfly,meanstatsperfly,controlstats,plottitle,'visible',visible);
  stathandles = PlotPerFrameStats2(stats_perframefeatures,statsperfly,meanstatsperfly,controlstats,plottitle,...
    'visible',visible,...
    'stattype','flymeans',...
    'weighttype','fracframesfly',...
    'plotstderr',true);
  drawnow;
  savename = sprintf('stats.png');
  savename = fullfile(figdir,savename);
  if exist(savename,'file'),
    delete(savename);
  end
  save2png(savename,stathandles.hfig);

  %stathandles = PlotPerFrameStatsExp(stats_perframefeatures,statsperexp,meanstatsperexp,controlstats,plottitle,'visible',visible);
  stathandles = PlotPerFrameStats2(stats_perframefeatures,statsperexp,meanstatsperexp,controlstats,plottitle,...
    'visible',visible,...
    'stattype','expmeans',...
    'plotstderr',true,'hfig',2);
    
  drawnow;
  savename = sprintf('statsexp.png');
  savename = fullfile(figdir,savename);
  if exist(savename,'file'),
    delete(savename);
  end
  save2png(savename,stathandles.hfig);
  
end

%% plot histograms

histfields = cell(1,numel(hist_perframefeatures));
histids = cell(1,numel(hist_perframefeatures));
for i = 1:numel(hist_perframefeatures),
  histfields{i} = hist_perframefeatures(i).field;
  if ismember(histfields{i},special_cases),
    frameconditionparams = DecodeConditions(hist_perframefeatures(i).framecondition,frameconditiondict);
    m = regexp(frameconditionparams(1:2:end),'^[^_]+_(.+)_labels$','once','tokens');
    tmp = find(~cellfun(@isempty,m),1);
    histids{i} = [histfields{i},'_',m{tmp}{1}];
  else
    if strcmp(hist_perframefeatures(i).framecondition,'any'),
      histids{i} = histfields{i};
    else
      histids{i} = [histfields{i},'_',hist_perframefeatures(i).framecondition];
    end
  end
end

% histfields = {hist_perframefeatures.field};
% histids = histfields;
% idxdur = find(ismember(histids,special_cases));
% for i = idxdur(:)',
%   frameconditionparams = DecodeConditions(hist_perframefeatures(i).framecondition,frameconditiondict);
%   m = regexp(frameconditionparams(1:2:end),'^[^_]+_(.+)_labels$','once','tokens');
%   tmp = find(~cellfun(@isempty,m),1);
%   histids{i} = [histids{i},'_',m{tmp}{1}];
% end

[histids,tmp,histidx] = unique(histids);
histfields = histfields(tmp);


for i = 1:numel(histfields),
  
  field = histfields{i};
  id = histids{i};
  idxcurr = find(histidx == i);
  if strcmp(field,'duration'),
    binfn = id;
  else
    binfn = field;
  end

  if ~isempty(controlhist),
    handles_control = PlotPerFrameHists3(id,field,idxcurr,hist_perframefeatures,...
      controlhist.meanhistperfly,controlhist.histperfly,...
      bins.(binfn),hist_plot_params,plottitle,...
      'visible',visible,'linestyle',':','stdstyle','errorbar',...
      'stattype','flymeans','meanweighttype','nframesfly',...
      'stdweighttype','fracframesfly',...
      'plotstderr',true);
    hax = handles_control.hax;
  else
    hax = [];
  end
  
%   handles = PlotPerFrameHists2(id,field,idxcurr,hist_perframefeatures,...
%     meanhistperfly,histperfly,...
%     bins.(binfn),hist_plot_params,plottitle,...
%     'visible',visible,...
%     'hax',hax);
  
  handles = PlotPerFrameHists3(id,field,idxcurr,hist_perframefeatures,...
    meanhistperfly,histperfly,...
    bins.(binfn),hist_plot_params,plottitle,...
    'visible',visible,...
    'hax',hax,...
    'stattype','flymeans','meanweighttype','nframesfly',...
    'stdweighttype','fracframesfly',...
    'plotstderr',true);
  
  if ~isempty(controlhist) && ~isempty(handles.hleg)
    % fix legend
    s = get(handles.hleg,'String');
    legend([handles_control.htype(1),handles.htype],[{'control'},s],'Parent',handles.hfig,'Interpreter','none');
  end

  
  drawnow;
  savename = sprintf('hist_%s.png',histfields{i});
  savename = fullfile(figdir,savename);
  if exist(savename,'file'),
    delete(savename);
  end
  save2png(savename,handles.hfig);
  
  if ~isempty(controlhist),
    handles_control = PlotPerFrameHist3(id,field,idxcurr,hist_perframefeatures,...
      controlhist.meanhistperexp,controlhist.histperexp,...
      bins.(binfn),hist_plot_params,plottitle,...
      'visible',visible,'linestyle',':','stdstyle','errorbar',...
      'stattype','expmeans',...
      'plotstderr',true);
    hax = handles_control.hax;
  else
    hax = [];
  end
  
%   handles = PlotPerFrameHists2Exp(id,field,idxcurr,hist_perframefeatures,...
%     meanhistperexp,histperexp,...
%     bins.(binfn),hist_plot_params,plottitle,...
%     'visible',visible,...
%     'hax',hax)
  
  handles = PlotPerFrameHists3(id,field,idxcurr,hist_perframefeatures,...
    meanhistperexp,histperexp,...
    bins.(binfn),hist_plot_params,plottitle,...
    'visible',visible,...
    'hax',hax,...
    'stattype','expmeans',...
    'plotstderr',true);

  
  if ~isempty(controlhist) && ~isempty(handles.hleg),
    % fix legend
    s = get(handles.hleg,'String');
    legend([handles_control.htype(1),handles.htype],[{'control'},s],'Parent',handles.hfig,'Interpreter','none');
  end

  
  drawnow;
  savename = sprintf('hist_%s_exp.png',histfields{i});
  savename = fullfile(figdir,savename);
  if exist(savename,'file'),
    delete(savename);
  end
  save2png(savename,handles.hfig);
  
  
end

close all;