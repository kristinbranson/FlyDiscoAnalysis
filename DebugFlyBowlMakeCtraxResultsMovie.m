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
  rootdir = 'E:\Data\FlyBowl\CtraxTest20110407';
else
  settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
  %rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
  rootdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20110804/LearnCtraxParams/expdirs';
end

%% parameters

analysis_protocol = '20110804';
params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol};

%%

expdirs = dir(fullfile(rootdir,'*_*'));
expdirs = cellfun(@(s) fullfile(rootdir,s),{expdirs([expdirs.isdir]).name},'UniformOutput',false);

for i = 3:numel(expdirs),
  expdir = expdirs{i};
  try
    FlyBowlMakeCtraxResultsMovie(expdir,params{:});
  catch ME,
    getReport(ME)
  end
end