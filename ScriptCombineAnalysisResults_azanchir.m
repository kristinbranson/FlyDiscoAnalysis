%% set up paths

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
screen_type = 'non_olympiad_azanchir_housing_CS_20120225';
%screen_type = 'non_olympiad_azanchir_housing_CS_20120204';
%screen_type = 'non_olympiad_azanchir_mating_galit_CS_20120211';
%screen_type = 'non_olympiad_azanchir_nicotine_mathias_berlin_20120211';
%analysis_protocol = sprintf('20120220_%s',screen_type);
analysis_protocol = sprintf('20120228_%s',screen_type);
rootdatadir = fullfile('/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax',analysis_protocol,'results');

datalocparamsfilestr = 'dataloc_params.txt';
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);


%% use SAGE to get experiment list

SAGEParams = SetSAGEParamsAzanchir(screen_type);


exp_metadata = SAGEListBowlExperiments(SAGEParams{:},...
  'checkflags',false,'removemissingdata',false,'rootdir',rootdatadir);
expdirs = {exp_metadata.file_system_path};
nexpdirs = numel(expdirs);
expnames = cell(size(expdirs));
for i = 1:nexpdirs,
  [~,expnames{i}] = fileparts(expdirs{i});
end


%% copy files from sciserv

if ismember(screen_type,{'non_olympiad_azanchir_housing_CS_20120225'}),
  for i = 1:numel(exp_metadata),
    fsp = exp_metadata(i).file_system_path;
    sspath = fullfile('/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data',expnames{i});
    SymbolicCopyExperimentDirectory(sspath,rootdatadir);    
  end  
end
%% create protocol names

protocolcats = cell(1,nexpdirs);
for i = 1:nexpdirs,
  protocolcats{i} = sprintf('%s__%s__%s',exp_metadata(i).protocol,exp_metadata(i).handling_protocol,exp_metadata(i).rearing_protocol);
end
[unique_protocolcats,~,protocolidx] = unique(protocolcats);

% get shorthands for all of these protocols
protocols = cell(1,numel(unique_protocolcats));
for i = 1:numel(unique_protocolcats),
  protocols{i} = input(sprintf('Name for protocols %s: ',unique_protocolcats{i}),'s');
end
% protocols = {'isolate_CS',...
%   'group_CS',...
%   'isolate_CS_plus_female',...
%   'group_CS_plus_female'};

protocolnames = protocols(protocolidx);

%% use file system to get experiment list
% 
% tmp = dir(rootdatadir);
% expdirs = {};
% exp_metadata = [];
% expnames = {};
% for i = 1:numel(tmp),
%   if ~tmp(i).isdir,
%     continue;
%   end
%   if ismember(tmp(i).name,{'.','..'}),
%     continue;
%   end
%   tmpexpdir = fullfile(rootdatadir,tmp(i).name);
%   tmpmetadata = parseExpDir(tmpexpdir);
%   if ~isempty(tmpmetadata),
%     exp_metadata = structappend(exp_metadata,tmpmetadata);
%     expdirs{end+1} = tmpexpdir; %#ok<SAGROW>
%     expnames{end+1} = tmp(i).name; %#ok<SAGROW>
%   end
% end
% 
% nexpdirs = numel(expdirs);

%% get protocol names from config files

% protocolnames = cell(1,nexpdirs);
% for i = 1:nexpdirs,
%   [~,tmp] = fileparts(ls(fullfile(expdirs{i},'*EP00015_Azanchir*.txt')));
%   protocolnames{i} = strrep(tmp,'FlyBowlDataCaptureParams_EP00015_Azanchir_','');
% end

%% combine experiments into protocols

[protocols,~,protocolidx] = unique(protocolnames);

% combine data
for i = 1:numel(protocols),
  idx = protocolidx == i;
  exps = cellfun(@(s) ['FlyBowl_',s],expnames(idx),'UniformOutput',false);
  % pecularity with SAGE REST interface
  if numel(exps) == 1,
    exps = exps{1};
  end
  outdir = fullfile(rootdatadir,protocols{i});
  FlyBowlCombineExperiments(rootdatadir,outdir,'controldatadirstr','',...
    'experiment_name',exps,'checkflags',false,'removemissingdata',false,...
    'analysis_protocol',analysis_protocol);
end

%% combine experiments to get "control"

outdir = fullfile(rootdatadir,'allconditions');
exps = cellfun(@(s) ['FlyBowl_',s],expnames,'UniformOutput',false);
FlyBowlCombineExperiments(rootdatadir,outdir,'controldatadirstr','',...
  'experiment_name',exps,'checkflags',false,'removemissingdata',false,...
  'analysis_protocol',analysis_protocol);
controlstats = load(fullfile(outdir,dataloc_params.statsperframematfilestr));
  
%% conditions to compare

comparisons = ...
  {{'housing_group_size',{'singlehoused','grouphoused03','grouphoused05','grouphoused10','grouphoused20'}},...
  {'isolation_length',{'isolated0','isolated1','isolated2','isolated3','isolated4'}},...
  {'recovery_length',{'recovery0','recovery02','recovery06','recovery12','recovery24'}},...
  {'all',protocols}};
  
for comparisoni = 1:numel(comparisons),
  
  comparison_name = comparisons{comparisoni}{1};
  protocols_curr = comparisons{comparisoni}{2};
  nprotocols_curr = numel(protocols_curr);

%% plot all stats together

statsperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statsperframefeaturesfilestr);
stats_perframefeatures = ReadStatsPerFrameFeatures(statsperframefeaturesfile);
statframeconditionsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statframeconditionfilestr);
frameconditiondict = ReadParams(statframeconditionsfile);
statflyconditionsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statflyconditionfilestr);
flyconditiondict = ReadParams(statflyconditionsfile);
statsparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statsparamsfilestr);
stats_params = ReadParams(statsparamsfile);

histplotparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histplotparamsfilestr);
hist_plot_params = ReadParams(histplotparamsfile);
histperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframefeaturesfilestr);
hist_perframefeatures = ReadHistPerFrameFeatures(histperframefeaturesfile);
histperframebinsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframebinsfilestr);
load(histperframebinsfile,'bins');

%nprotocols = numel(protocols);
allmeanstatsperfly = cell(1,nprotocols_curr);
allstatsperfly = cell(1,nprotocols_curr);
for i = 1:nprotocols_curr,
  datacurr = load(fullfile(rootdatadir,protocols_curr{i},dataloc_params.statsperframematfilestr));
  allstatsperfly{i} = datacurr.statsperfly;
  allmeanstatsperfly{i} = datacurr.meanstatsperfly;
end

handles = PlotPerFrameStatsComparison(stats_perframefeatures,allstatsperfly,allmeanstatsperfly,controlstats,protocols_curr);
set(handles.hfig,'Visible','on');
filename = 'stats';
exts = {'pdf','jpeg'};
if strcmpi(comparison_name,'all'),
  outdir = fullfile(rootdatadir,'analysis_plots');
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

%% plot all histograms together

allmeanhistperfly = cell(1,nprotocols_curr);
allhistperfly = cell(1,nprotocols_curr);
for i = 1:nprotocols_curr,
  datacurr = load(fullfile(rootdatadir,protocols_curr{i},dataloc_params.histperframematfilestr));
  allhistperfly{i} = datacurr.histperfly;
  allmeanhistperfly{i} = datacurr.meanhistperfly;
end

%exts = {'jpeg','pdf'};
%outdir = fullfile(rootdatadir,'analysis_plots');

for typei = 1:numel(hist_perframefeatures),
  handles = PlotPerFrameHistsComparison(hist_perframefeatures(typei).field,...
    hist_perframefeatures(typei).flycondition,...
    hist_perframefeatures(typei).framecondition,...
    allmeanhistperfly,allhistperfly,...
    bins.(hist_perframefeatures(typei).field),...
    hist_plot_params,...
    protocols_curr,...
    'visible','on');
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
end
% 
% handles = PlotPerFrameHists(field,hist_perframefeatures,...
%     meanhistperfly,histperfly,...
%     bins.(field),hist_plot_params,plottitle,...
%     'visible',visible,...
%     'hax',hax);

%% end loop over comparisons

end