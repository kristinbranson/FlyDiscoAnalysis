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
linename = 'pBDPGAL4U';
requiredfiles = {'histperframematfilestr','statsperframematfilestr'};
fn = 'velmag_flyany_frameany';

%% get all experiments that satisfy input conditions; currently parsing
% directory structure
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
subreadfiles = cell(size(requiredfiles));
for i = 1:numel(requiredfiles),
  subreadfiles{i} = dataloc_params.(requiredfiles{i});
end
subreadfiles = unique(subreadfiles);

exp_params = {'rootdir',rootdir,'subreadfiles',subreadfiles,'linename',linename};

[expdirs,expdir_reads,~,experiments,~,~] = ...
  getExperimentDirs(exp_params{:});
nexpdirs = numel(expdirs);

if nexpdirs == 0,
  error('No experiments selected');
end


%%
hfig = CompareHistograms(expdir_reads,fn);