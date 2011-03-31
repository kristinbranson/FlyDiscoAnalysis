%% set up paths

if ispc,
  addpath E:\Code\JCtrax\misc;
  addpath E:\Code\JCtrax\filehandling;
  settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
  %expdir = 'E:\Data\FlyBowl\bowl_data\pBDPGAL4U_TrpA_Rig1Plate10BowlB_20110217T145323';
  expdir = 'E:\Data\FlyBowl\bowl_data\GMR_14B02_AE_01_TrpA_Rig2Plate14BowlD_20110204T094327';
  addpath E:\Code\Ctrax\trunk\matlab\netlab;
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
  addpath /groups/branson/home/bransonk/tracking/code/Ctrax/matlab/netlab;
  settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
  expdir = '../fly_bowl_sciserv/bowl_data/GMR_13B10_AE_01_TrpA_Rig2Plate14BowlC_20110204T091600';
end

%% parameters

analysis_protocol = '20110222';

params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol};

%% 
FlyBowlExtraDiagnostics(expdir,params{:});
