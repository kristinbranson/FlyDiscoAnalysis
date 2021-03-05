%DebugFlyDiscoRegisterTrx.m 4/15/2016
% 
% /groups/flyprojects/home/leea30/git/fba.build/bubble/current/run_FlyDiscoRegisterTrx.sh /groups/branson/bransonlab/projects/olympiad/MCR/v717 /tier2/branson/fly_bubble/01_tracked/social_GMR_22D03_AE_01_dTrpA1_RigA_20160325T134540 
% settingsdir = ' /groups/flyprojects/home/olympiad/bubble_bin/FlyBubbleAnalysis/settings/ analysis_protocol current_bubble'
%%
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings';
%rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
rootdir = '/tier2/branson/FlyBubble_testing/temp/testing_registration';
expdir = fullfile(rootdir,'grooming_GMR_60E02_AE_01_CsChr_RigA_20160519T164300');

%% parameters

% analysis_protocol = 'current_bubble';
analysis_protocol = 'current_non_olympiad_rubin_grooming';
params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol};

%% run once

FlyDiscoRegisterTrx(expdir,params{:});