% test FlyBowlRegisterTrx

%% set up path
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;

%% get experiments

analysis_protocol = '20110202';

[expdirs,expdir_reads,expdir_writes,experiments,rootreaddir,rootwritedir] = ...
  getExperimentDirs('protocol',analysis_protocol,...
  'subreadfiles',{'movie.ufmf.ann','ctrax_results.mat'});

%% do the registration

for i = 1:numel(expdirs),
  fprintf('Registering experiment %s\n',expdirs{i});
  FlyBowlRegisterTrx(expdir_reads{i},'analysis_protocol',analysis_protocol);
end