%% set up paths

bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
end

%% data locations

if ispc,
  settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
  rootdir = 'E:\Data\FlyBowl\bowl_data';
else
  settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
  %rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
  rootdir = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/CtraxTest20110407';
end

%% parameters

analysis_protocol = '20110407';
params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol,...
  'visible','on'};

%%

expdirs = dir(fullfile(rootdir,'*_*'));
expdirs = cellfun(@(s) fullfile(rootdir,s),{expdirs([expdirs.isdir]).name},'UniformOutput',false);

for i = 1:numel(expdirs),
  expdir = expdirs{i};
  FlyBowlMakeCtraxResultsMovie(expdir,params{:});
end