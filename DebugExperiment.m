% DebugExperiment

%% paths

if ispc
  
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
end

%% locations of data

if ispc,
  
  rootreaddir = 'E:\Data\FlyBowl\polycarbonate_scratched';
  rootwritedir = 'E:\Data\FlyBowl\polycarbonate_scratched';
  expdir = 'DL-wildtype_TrpA_Rig1Plate01BowlC_20101112T152624';
  addpath E:\Code\JCtrax\misc;
  addpath E:\Code\JCtrax\filehandling;
  nr = 1024;
  nc = 1024;
  ncolors = 1;
  nframes = 9053;
  vidinfo = struct('nr',nr,'nc',nc,'ncolors',ncolors,'nframes',nframes);
else
  protocol = 'CtraxTest20110111';
  [expdirs,expdir_reads,expdir_writes,experiments,rootreaddir,rootwritedir] = ...
    getExperimentDirs('protocol',protocol);
  expdir = expdirs{1};
  [readframe,nframes,fid,vidinfo] = get_readframe_fcn(fullfile(expdir_reads{1},'movie.ufmf'));
  fclose(fid);
end

%%

obj = Experiment('rootreaddir',rootreaddir,'rootwritedir',rootwritedir);
obj.AddExpDir(expdir,vidinfo);