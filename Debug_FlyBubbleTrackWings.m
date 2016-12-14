% Debug wing tracking and choose orientations based on wing tracking
% 12/13/2016

addpath /groups/branson/home/robiea/Code_versioned/Jdetect_github/Jdetect/filehandling/;
addpath /groups/branson/home/robiea/Code_versioned/Jdetect_github/Jdetect/misc;
addpath /groups/branson/home/robiea/Code_versioned/JCtrax/simple_wing;
settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings';
%rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
rootdir = '/tier2/branson/FlyBubble_testing/temp/testing_registration';
% expdir = fullfile(rootdir,'grooming_GMR_60E02_AE_01_CsChr_RigA_20160519T164300');
expdir = '/groups/branson/home/robiea/Projects_data/temp/cx_GMR_OL0077B_CsChr_RigB_20150901T133146';

%% parameters

analysis_protocol = 'current_bubble';
% analysis_protocol = 'current_non_olympiad_rubin_grooming';
params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol};

%% run once



[trx,perframedata,info,wtunits,trackdata] = FlyBubbleTrackWings(expdir,params{:});