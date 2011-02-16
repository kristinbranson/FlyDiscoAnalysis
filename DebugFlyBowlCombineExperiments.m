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
analysis_protocol = '20110211';
linename = 'pBDPGAL4U';
outdir = fullfile('/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/line_data/pBDPGAL4U',datestr(now,30));

%%

FlyBowlCombineExperiments(rootdir,outdir,...
  'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol,...
  'linename',linename);
