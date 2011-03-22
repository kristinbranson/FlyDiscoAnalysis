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
  rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';  
end
analysis_protocol = '20110222';
expdir = fullfile(rootdir,'GMR_14B02_AE_01_TrpA_Rig2Plate14BowlD_20110204T094327');

%% run

FlyBowlMakeCtraxResultsMovie(expdir,...
  'analysis_protocol',analysis_protocol,'settingsdir',settingsdir);