% combine data from multiple experiments
function FlyBowlCombineExperiments(rootdir,outdir,varargin)

% parse arguments
[settingsdir,analysis_protocol,datalocparamsfilestr,...
  requiredfiles,subreadfiles,visible,leftovers] = ...
  myparse_nocheck(varargin,...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'analysis_protocol','current',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'requiredfiles',{'histperframematfilestr','statsperframematfilestr'},...
  'subreadfiles',{},...
  'visible','on');

%% get all experiments that satisfy input conditions; currently parsing
% directory structure
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
for i = 1:numel(requiredfiles),
  subreadfiles{end+1} = dataloc_params.(requiredfiles{i}); %#ok<AGROW>
end
subreadfiles = unique(subreadfiles);

exp_params = [{'rootdir',rootdir,'subreadfiles',subreadfiles},leftovers];

[expdirs,expdir_reads,~,experiments,~,~] = ...
  getExperimentDirs(exp_params{:});
nexpdirs = numel(expdirs);

if nexpdirs == 0,
  error('No experiments selected');
end

%% load per-fly data

statsperfly = [];
histperfly = [];

for i = 1:nexpdirs,
  
  % load the data for the current experiment
  %i=i+1
  expdir = expdir_reads{i};
  statsmatfile = fullfile(expdir,dataloc_params.statsperframematfilestr);
  statsdata = load(statsmatfile);
  histmatfile = fullfile(expdir,dataloc_params.histperframematfilestr);
  histdata = load(histmatfile);

  % merge with data read in so far
  statsperfly = MergeStatsPerFly(statsperfly,statsdata.statsperfly,expdir);
  %numel(statsperfly.a_mm_flyany_frameany.Z)
  histperfly = MergeHistPerFly(histperfly,histdata.histperfly,expdir);

end

%% Combine per-fly data as in one experiment

statsperexp = struct;
histperexp = struct;

fns = fieldnames(statsperfly);
for i = 1:numel(fns),
  fn = fns{i};

  % compute per-experiment hist
  statsperexp.(fn) = CombinePerFrameStats(statsperfly.(fn));
  
end

fns = fieldnames(histperfly);
for i = 1:numel(fns),
  fn = fns{i};

  % compute per-experiment hist
  histperexp.(fn) = CombinePerFrameHists(histperfly.(fn));
  
end

%% create output directory
if ~exist(outdir,'file'),
  [success,msg] = mkdir(outdir);
  if ~success,
    error('Could not create directory %s for output:\n%s',outdir,msg);
  end
end

%% save 

statstxtsavename = fullfile(outdir,dataloc_params.statsperframetxtfilestr);
statsmatsavename = fullfile(outdir,dataloc_params.statsperframematfilestr);

statsdata.statsperfly = statsperfly;
statsdata.statsperexp = statsperexp;
statsdata.exp_params = exp_params;
statsdata.expdirs = expdir_reads;
statsdata.experiments = experiments;

save(statsmatsavename,'-struct','statsdata');
% save to text file
fns = fieldnames(statsperfly);
for i = 1:numel(fns),
  fn = fns{i};
  SavePerFrameStatsTxtFile(statstxtsavename,fn,statsperfly.(fn),statsperexp.(fn));
end

histtxtsavename = fullfile(outdir,dataloc_params.histperframetxtfilestr);
histmatsavename = fullfile(outdir,dataloc_params.histperframematfilestr);

histdata.histperfly = histperfly;
histdata.histperexp = histperexp;
histdata.exp_params = exp_params;
histdata.expdirs = expdir_reads;
histdata.experiments = experiments;

save(histmatsavename,'-struct','histdata');
% save to text file
fns = fieldnames(histperfly);
for i = 1:numel(fns),
  fn = fns{i};
  SavePerFrameHistTxtFile(histtxtsavename,fn,histperfly.(fn),histperexp.(fn));
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
fprintf(fid,',%s',expdir_reads{:});
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
[~,typename] = fileparts(outdir);
histperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframefeaturesfilestr);
hist_perframefeatures = ReadHistPerFrameFeatures(histperframefeaturesfile);
histperframebinsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframebinsfilestr);
load(histperframebinsfile,'bins');

%% plot stuff
% TODO: plot stats of controls for the previous X weeks behind these

hist_fields = unique({hist_perframefeatures.field});
for i = 1:numel(hist_fields),
  
  field = hist_fields{i};

  hfig = PlotPerFrameHists(field,hist_perframefeatures,...
    histperexp,histperfly,...
    bins.(field),hist_plot_params,typename,...
    'visible',visible);
  drawnow;
  savename = fullfile(figdir,sprintf('hist_%s.png',hist_fields{i}));
  if exist(savename,'file'),
    delete(savename);
  end
  save2png(savename,hfig);
  
end

close all;

