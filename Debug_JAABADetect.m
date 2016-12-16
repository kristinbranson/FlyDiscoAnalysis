%Debug_JAABADetect.m 12/16/2016

%% setup path

addpath /groups/branson/home/robiea/Code_versioned/Jdetect_github/Jdetect/filehandling/;
addpath /groups/branson/home/robiea/Code_versioned/Jdetect_github/Jdetect/misc;
addpath /groups/branson/home/robiea/Code_versioned/Jdetect_github/Jdetect/perframe/;

settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings';

%% experiment directories
expdir = '/groups/branson/home/robiea/Projects_data/temp/grooming2_Line1_Unknown_RigA_20161213T224309';

%% parameters
analysis_protocol = '20161215_flybubble_seedsa';
% analysis_protocol = 'current_non_olympiad_rubin_grooming';
params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol};
%%

datalocparamsfilestr = 'dataloc_params.txt';
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
jaabaclassifierparamsfilestrs = fullfile(settingsdir,analysis_protocol,dataloc_params.jaabaclassifierparamsfilestrs);

jabfiles = textread(jaabaclassifierparamsfilestrs,'%s');

pwdprev = pwd;
jaabadir = fileparts(which('JAABADetect'));
cd(jaabadir);

JAABADetect(expdir,'jabfiles',jabfiles,'forcecompute',false);

cd(pwdprev);