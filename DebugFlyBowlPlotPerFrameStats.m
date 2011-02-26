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
  analysis_protocol = '20110202_pc';
  dataloc_params = ReadParams(fullfile(settingsdir,analysis_protocol,'dataloc_params.txt'));
  expdir = 'pBDPGAL4U_TrpA_Rig1Plate10BowlA_20110202T105734';
  expdir_read = fullfile(dataloc_params.rootreaddir,expdir);
  
else
  settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
  protocol = '';
  analysis_protocol = '20110211';
  [~,expdirs_ready] = ...
    getExperimentDirs('protocol',analysis_protocol,'subreadfiles',{'hist_perframe.mat'});
  [~,expdirs_done] = ...
    getExperimentDirs('protocol',analysis_protocol,'subreadfiles',{'analysis_plots/hist_velmag.png'});
  expdirs = setdiff(expdirs_ready,expdirs_done);
  expdir = expdirs{1};
end

%% run

FlyBowlPlotPerFrameStats(expdir,'visible','on');