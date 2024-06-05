%DebugFlyDiscoDetectIndicatorLedOnOff.m 4/25/2016
% 
% /groups/flyprojects/home/leea30/git/fba.build/bubble/current/run_FlyDiscoRegisterTrx.sh /groups/branson/bransonlab/projects/olympiad/MCR/v717 /tier2/branson/fly_bubble/01_tracked/social_GMR_22D03_AE_01_dTrpA1_RigA_20160325T134540 
% settingsdir = ' /groups/flyprojects/home/olympiad/bubble_bin/FlyBubbleAnalysis/settings/ analysis_protocol current_bubble'
%
modpath
settingsdir = '/home/robiea@hhmi.org/Code_versioned/FlyDiscoAnalysis/settings-internal';
%rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
rootdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240507_testingVNC3';
expdir = fullfile(rootdir,'VNC3_YNA_K_162984_RigD_20240507T182533');

%% parameters

% analysis_protocol = 'current_bubble';
% analysis_protocol = '20160415_flybubble_non_optogenetic';
analysis_protocol = '20240507_flybubble_LED_VNC3';
params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol};

%% run once

FlyDiscoDetectIndicatorLedOnOff(expdir,params{:});