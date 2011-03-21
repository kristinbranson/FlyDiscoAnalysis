if ispc,
  addpath E:\Code\JCtrax\misc;
  addpath E:\Code\JCtrax\filehandling;
  settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
  rootdir = 'E:\Data\FlyBowl\bowl_data';
  expdir = 'E:\Data\FlyBowl\bowl_data\GMR_14B02_AE_01_TrpA_Rig2Plate14BowlD_20110204T094327';
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
  rootdir = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/CtraxTest20110315';
  expdir = fullfile(rootdir,'GMR_71C10_AE_01_TrpA_Rig1Plate10BowlA_20110311T152418');
end

%% parameters

analysis_protocol = '20110222';
params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol};

%% run once

FlyBowlRegisterTrx(expdir,params{:});

%% run a bunch of times

expdirs = dir(fullfile(rootdir,'*_*'));
expdirs = {expdirs.name};
for i = 1:numel(expdirs),
  FlyBowlRegisterTrx(fullfile(rootdir,expdirs{i}),params{:});
  drawnow;
end