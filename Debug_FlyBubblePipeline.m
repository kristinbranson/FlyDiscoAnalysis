%Debug_FlyBubblePipeline.m

addpath /groups/branson/home/robiea/Code_versioned/Jdetect_github/Jdetect/filehandling/;
addpath /groups/branson/home/robiea/Code_versioned/Jdetect_github/Jdetect/misc;
addpath /groups/branson/home/robiea/Code_versioned/JCtrax/simplewing;
addpath //groups/branson/home/bransonk/tracking/code/lds/hmm;
addpath /groups/branson/home/robiea/Code_versioned/flySpaceTimeFeatures;
addpath /groups/branson/home/robiea/Code_versioned/Jdetect_github/Jdetect/perframe;
settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings';

%% parameters
analysis_protocol = '20161215_flybubble_seedsa';
% analysis_protocol = 'current_non_olympiad_rubin_grooming';
params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol};

%% experiment directory
expdir = '/groups/branson/home/robiea/Projects_data/temp/testing_pipelinescript';

%% run

[success,msgs,stage] = FlyBubblePipeline(expdir,params{:})