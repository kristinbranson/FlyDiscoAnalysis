%Debug_preparePerFrameFtrs 12/15/2016

%% setup path
% code from  https://github.com/mkabra/flySpaceTimeFeatures
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis;
% addpath /groups/branson/home/robiea/Code_versioned/flySpaceTimeFeatures;
% addpath /groups/branson/home/robiea/Code_versioned/flySpaceTimeFeatures/classifier;
% addpath /groups/branson/home/robiea/Code_versioned/flySpaceTimeFeatures/figurecode;
% addpath /groups/branson/home/robiea/Code_versioned/flySpaceTimeFeatures/deepmatching;
cd /groups/branson/home/robiea/Code_versioned/flySpaceTimeFeatures

%% set experiments

expdir = '/groups/branson/home/robiea/Projects_data/temp/cx_GMR_OL0077B_CsChr_RigB_20150901T133146';
%% set parameters
settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings';
analysis_protocol = 'current_bubble';
datalocparamsfilestr = 'dataloc_params.txt';
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
trxfilestr = fullfile(expdir,dataloc_params.trxfilestr);
movfilestr = fullfile(expdir,dataloc_params.moviefilestr);


%% run code
% preparePerFrameFtrs(pathToMovie,pathToTrxfile, false, false);
%flies_hoghof_hs_notaligned -- flow computed using Horn-Schunck. The current frame and next frame are not aligned in any way.

preparePerFrameFtrs(movfilestr,trxfilestr,false,false);
