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
  rootdir = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/CtraxTest20110315';
end
analysis_protocol = '20110222';
expdir = fullfile(rootdir,'GMR_71C10_AE_01_TrpA_Rig1Plate10BowlA_20110311T152418');

params = {'analysis_protocol',analysis_protocol,'settingsdir',settingsdir};

%% run

FlyBowlMakeCtraxResultsMovie(expdir,params{:});