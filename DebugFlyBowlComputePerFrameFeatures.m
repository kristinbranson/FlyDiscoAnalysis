%% set up paths

if ispc,
  addpath E:\Code\JCtrax\misc;
  addpath E:\Code\JCtrax\filehandling;
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
end

%% data locations

if ispc,
  settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
  rootdir = 'E:\Data\FlyBowl\bowl_data';
else
  settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
  %rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
  rootdir = '/groups/branson/bransonlab/tracking_data/olympiad/HackHitData';
end

%% parameters

analysis_protocol = 'current';
params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol};

%% 
rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
experiment_name = 'EXT_CantonS_1220002_None_Rig1Plate15BowlA_20120519T124001';
expdir = fullfile(rootdatadir,experiment_name);
SymbolicCopyExperimentDirectory(expdir,rootdir);
newexpdir = fullfile(rootdir,experiment_name);
FlyBowlComputePerFrameFeatures(newexpdir,params{:});

%%

expdirs = dir(fullfile(rootdir,'*_*'));
expdirs = cellfun(@(s) fullfile(rootdir,s),{expdirs([expdirs.isdir]).name},'UniformOutput',false);

for i = 1:numel(expdirs),
  FlyBowlComputePerFrameFeatures(expdirs{i},params{:});
end
