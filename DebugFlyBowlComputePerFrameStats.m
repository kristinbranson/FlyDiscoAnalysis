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
  protocol = 'CtraxTest20110202';
  analysis_protocol = '20110211';
  [expdirs,expdir_reads,expdir_writes,experiments,rootreaddir,rootwritedir] = ...
    getExperimentDirs('protocol',protocol,'subreadfiles',{'registered_trx.mat'});
end
% [readframe,nframes,fid,vidinfo] = get_readframe_fcn(fullfile(expdir_read,'movie.ufmf'));
% fclose(fid);

%% run

for expdiri = randperm(numel(expdirs)),
  expdir = expdir_reads{expdiri};
  FlyBowlComputePerFrameStats(expdir,'analysis_protocol',analysis_protocol);
end