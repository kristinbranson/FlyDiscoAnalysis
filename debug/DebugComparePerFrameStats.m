% DebugFlyBowlCombineExperiments

% set up path
if ispc,
  addpath E:\Code\JCtrax\misc;
  addpath E:\Code\JCtrax\filehandling;
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
end
 
%% locations of data

rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
datalocparamsfilestr = 'dataloc_params.txt';
analysis_protocol = '20110211';
linename1 = 'pBDPGAL4U';
linename2 = 'GMR_15C06_AE_01';
requiredfiles = {'statsperframematfilestr'};
field = 'velmag';

%% get all experiments that satisfy input conditions; currently parsing
% directory structure
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
subreadfiles = cell(size(requiredfiles));
for i = 1:numel(requiredfiles),
  subreadfiles{i} = dataloc_params.(requiredfiles{i});
end
subreadfiles = unique(subreadfiles);

%% experiments for line 1
exp_params = {'rootdir',rootdir,'subreadfiles',subreadfiles,'linename',linename1};
[~,expdirs1,~,~,~,~] = ...
  getExperimentDirs(exp_params{:});
nexpdirs1 = numel(expdirs1);

if nexpdirs1 == 0,
  error('No experiments selected');
end

%% experiments for line 2
exp_params = {'rootdir',rootdir,'subreadfiles',subreadfiles,'linename',linename2};
[~,expdirs2,~,~,~,~] = ...
  getExperimentDirs(exp_params{:});
nexpdirs2 = numel(expdirs2);

if nexpdirs2 == 0,
  error('No experiments selected');
end


%%
colors = cat(1,repmat([0,0,0],[nexpdirs1,1]),...
  repmat([.7,0,0],[nexpdirs2,1]));

hfig = ComparePerFrameStats([expdirs1,expdirs2],field,'colors',colors,'sortby','none');