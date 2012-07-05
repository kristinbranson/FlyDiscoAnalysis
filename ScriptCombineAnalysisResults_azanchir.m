%% set up paths

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
screen_type = 'non_olympiad_azanchir_housing_CS_20120225';
%screen_type = 'non_olympiad_azanchir_housing_CS_20120204';
screen_type = 'non_olympiad_azanchir_mating_galit_CS_20120211';
%screen_type = 'non_olympiad_azanchir_nicotine_mathias_berlin_20120211';
analysis_protocol = sprintf('20120220_%s',screen_type);
%analysis_protocol = sprintf('20120228_%s',screen_type);
rootdatadir = fullfile('/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax',analysis_protocol,'results');

%analysis_protocol = '20120330';

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

%% copy over behavior labels

tmpinputdir = '/groups/branson/home/kabram/clusterScripts/matlaboutWalk';
for i = 1:numel(expdirs),
  expdir = expdirs{i};
  [~,experiment_name] = fileparts(expdir);
  inexpdir = fullfile(tmpinputdir,experiment_name);
  if ~exist(inexpdir,'dir'),
    warning('No directory %s containing behavior labels',inexpdir);
    continue;
  end
  tmp = dir(fullfile(inexpdir,'*labels*.mat'));
  fprintf('%d label files found for %s\n',numel(tmp),experiment_name);
  for j = 1:numel(tmp),
    cmd = sprintf('cp %s %s',fullfile(inexpdir,tmp(j).name),expdir);
    unix(cmd);
  end
  tmp = dir(fullfile(inexpdir,'*scores*.mat'));
  fprintf('%d score files found for %s\n',numel(tmp),experiment_name);
  for j = 1:numel(tmp),
    cmd = sprintf('cp %s %s',fullfile(inexpdir,tmp(j).name),expdir);
    unix(cmd);
  end
end

%% compute new per-frame stats

for i = 1:numel(expdirs),
  try
  FlyBowlComputePerFrameStats2(expdirs{i},'analysis_protocol',analysis_protocol);
  FlyBowlPlotPerFrameStats2(expdirs{i},'analysis_protocol',analysis_protocol,'controldatadirstr','');
  catch ME
    getReport(ME)
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

% Name for protocols EP_flybowl_v011p3.xls__HP_flybowl_v007p1.xls__RP_Olympiad_v010p0.xls: mated
% Name for protocols EP_flybowl_v011p3.xls__HP_flybowl_v007p2.xls__RP_Olympiad_v010p0.xls: rejected
% Name for protocols EP_flybowl_v011p3.xls__HP_flybowl_v007p3.xls__RP_Olympiad_v010p0.xls: naive
% Name for protocols EP_flybowl_v011p3.xls__HP_flybowl_v007p6.xls__RP_Olympiad_v010p0.xls: 
% Name for protocols EP_flybowl_v011p3.xls__HP_flybowl_v007p7.xls__RP_Olympiad_v010p0.xls: 

% Name for protocols EP_flybowl_v011p1.xls__HP_flybowl_v007p8.xls__RP_Olympiad_v010p0.xls: grouphoused20
% Name for protocols EP_flybowl_v011p1.xls__HP_flybowl_v007p9.xls__RP_Olympiad_v010p0.xls: grouphoused10
% Name for protocols EP_flybowl_v011p1.xls__HP_flybowl_v008p0.xls__RP_Olympiad_v010p0.xls: grouphoused05
% Name for protocols EP_flybowl_v011p1.xls__HP_flybowl_v008p1.xls__RP_Olympiad_v010p0.xls: grouphoused03
% Name for protocols EP_flybowl_v011p1.xls__HP_flybowl_v008p2.xls__RP_Olympiad_v010p0.xls: singlehoused
% Name for protocols EP_flybowl_v011p1.xls__HP_flybowl_v008p3.xls__RP_Olympiad_v010p0.xls: isolated0
% Name for protocols EP_flybowl_v011p1.xls__HP_flybowl_v008p4.xls__RP_Olympiad_v010p0.xls: isolated4
% Name for protocols EP_flybowl_v011p1.xls__HP_flybowl_v008p5.xls__RP_Olympiad_v010p0.xls: isolated3
% Name for protocols EP_flybowl_v011p1.xls__HP_flybowl_v008p6.xls__RP_Olympiad_v010p0.xls: isolated2
% Name for protocols EP_flybowl_v011p1.xls__HP_flybowl_v008p7.xls__RP_Olympiad_v010p0.xls: isolated1
% Name for protocols EP_flybowl_v011p1.xls__HP_flybowl_v008p8.xls__RP_Olympiad_v010p0.xls: recovery00
% Name for protocols EP_flybowl_v011p1.xls__HP_flybowl_v008p9.xls__RP_Olympiad_v010p0.xls: recovery24
% Name for protocols EP_flybowl_v011p1.xls__HP_flybowl_v009p0.xls__RP_Olympiad_v010p0.xls: recovery12
% Name for protocols EP_flybowl_v011p1.xls__HP_flybowl_v009p1.xls__RP_Olympiad_v010p0.xls: recovery06
% Name for protocols EP_flybowl_v011p1.xls__HP_flybowl_v009p2.xls__RP_Olympiad_v010p0.xls: recovery02

protocolnames = protocols(protocolidx);
expignore = strcmp(protocolnames,'');
protocols(strcmp(protocols,'')) = [];

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

%% combine experiments to get "control"

outdir = fullfile(rootdatadir,'allconditions');
exps = cellfun(@(s) ['FlyBowl_',s],expnames(~expignore),'UniformOutput',false);
FlyBowlCombineExperiments2(rootdatadir,outdir,'controldatadirstr','',...
  'experiment_name',exps,'checkflags',false,'removemissingdata',false,...
  'analysis_protocol',analysis_protocol,...
  'visible','off');
controlstats = load(fullfile(outdir,dataloc_params.statsperframematfilestr));

%% combine experiments into protocols

% combine data
for i = 1:numel(protocols),
  if strcmp(protocols,''),
    continue;
  end
  idx = strcmp(protocolnames,protocols{i});
  %idx = protocolidx == i;
  exps = cellfun(@(s) ['FlyBowl_',s],expnames(idx),'UniformOutput',false);
  % pecularity with SAGE REST interface
  if numel(exps) == 1,
    exps = exps{1};
  end
  outdir = fullfile(rootdatadir,protocols{i});
  FlyBowlCombineExperiments2(rootdatadir,outdir,'controldatadirstr','',...
    'experiment_name',exps,'checkflags',false,'removemissingdata',false,...
    'analysis_protocol',analysis_protocol,...
    'controlstats',controlstats,...
    'visible','off');
end


  
%% conditions to compare
% 
% comparisons = ...
%   {{'housing_group_size',{'singlehoused','grouphoused03','grouphoused05','grouphoused10','grouphoused20'}},...
%   {'isolation_length',{'isolated0','isolated1','isolated2','isolated3','isolated4'}},...
%   {'recovery_length',{'recovery00','recovery02','recovery06','recovery12','recovery24'}},...
%   {'all',protocols}};

comparisons = ...
  {{'all',protocols}};

%%

for comparisoni = 1:numel(comparisons),

  fprintf('Plotting data for comparison %s...\n',comparisons{comparisoni}{1});
  FlyBowlPlotComparisons(rootdatadir,comparisons{comparisoni}{2},comparisons{comparisoni}{1},...
    'analysis_protocol',analysis_protocol,'settingsdir',settingsdir,...
    'controlstats',controlstats,'fly_plotstderr',true,...
    'exp_plotstderr',true,...
    'meanweighttype','nframesfly',...
    'stdweighttype','fracframesfly');
end

