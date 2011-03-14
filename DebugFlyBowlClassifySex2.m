%% set up path

if ispc,
  addpath E:\Code\JCtrax\misc;
  addpath E:\Code\JCtrax\filehandling;
  addpath E:\Code\hmm;
  settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
  rootdir = 'E:\Data\FlyBowl\bowl_data';
  expdir = 'E:\Data\FlyBowl\bowl_data\GMR_14B02_AE_01_TrpA_Rig2Plate14BowlD_20110204T094327';
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/lds/hmm;
  settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
  rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
  expdir = fullfile(rootdir,'GMR_12E07_AE_01_TrpA_Rig2Plate14BowlA_20110209T134320');
end

%% parameters

analysis_protocol = '20110222';
params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol};

%% 

FlyBowlClassifySex2(expdir,params{:});