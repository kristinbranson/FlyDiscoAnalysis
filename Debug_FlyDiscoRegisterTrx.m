%DebugFlyDiscoRegisterTrx.m 4/15/2016
% 
% /groups/flyprojects/home/leea30/git/fba.build/bubble/current/run_FlyDiscoRegisterTrx.sh /groups/branson/bransonlab/projects/olympiad/MCR/v717 /tier2/branson/fly_bubble/01_tracked/social_GMR_22D03_AE_01_dTrpA1_RigA_20160325T134540 
% settingsdir = ' /groups/flyprojects/home/olympiad/bubble_bin/FlyBubbleAnalysis/settings/ analysis_protocol current_bubble'
%% set path
% addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
% addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
% settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings';
% %rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
% rootdir = '/tier2/branson/FlyBubble_testing/temp/testing_registration';
modpath
%%
settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/settings';
% expdir = fullfile(rootdir,'grooming_GMR_60E02_AE_01_CsChr_RigA_20160519T164300');

%% parameters

% analysis_protocol = 'current_bubble';
% analysis_protocol = 'current_non_olympiad_rubin_grooming';
% analysis_protocol = '20210531_flybubble_LED_ARtestingperframe';
analysis_protocol = '20211014_flybubbleRed_noChr';
params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol};

%% test expdir
% expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210802_testACI_movielength/VNC_JRC_SS65710_RigB_20210513T142817';
% expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20211014_testing_FlyBubbleRed/socialCsChr_GMR_72C11_AE_01_CsChrimson_RigD_20191114T172654';
expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20211014_testing_FlyBubbleRed/nochr_TrpApBDP_Unknown_RigB_20201216T160731';

%% run once

FlyDiscoRegisterTrx(expdir,params{:});