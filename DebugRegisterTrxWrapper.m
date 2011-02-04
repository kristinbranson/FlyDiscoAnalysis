% DebugRegisterTrxWrapper

%% set up path

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;

%% data locations

protocol = '20110111';
[expdirs,expdir_reads,expdir_writes,experiments,rootreaddir,rootwritedir] = ...
  getExperimentDirs('protocol',['CtraxTest',protocol]);

%% register all experiments

for i = 1:numel(expdirs),
  RegisterTrxWrapper(expdirs{i},protocol);
end