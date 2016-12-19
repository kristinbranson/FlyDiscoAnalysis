% Debug_FlyBubbleAutomaticChecks_Complete 12/19/2016


addpath /groups/branson/home/robiea/Code_versioned/Jdetect_github/Jdetect/filehandling/;
addpath /groups/branson/home/robiea/Code_versioned/Jdetect_github/Jdetect/misc;
addpath /groups/branson/home/robiea/Code_versioned/Jdetect_github/Jdetect/perframe/;

settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings';

%% parameters
analysis_protocol = '20161215_flybubble_seedsa';
% analysis_protocol = 'current_non_olympiad_rubin_grooming';
params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol};

%% experiment directory
expdir = '/groups/branson/home/robiea/Projects_data/temp/grooming2_Line1_Unknown_RigA_20161213T224309';

%% run

FlyBubbleAutomaticChecks_Complete(expdir,params{:});