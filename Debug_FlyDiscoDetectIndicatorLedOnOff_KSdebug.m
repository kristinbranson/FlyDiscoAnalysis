%DebugFlyDiscoDetectIndicatorLedOnOff.m 4/25/2016
% 
% /groups/flyprojects/home/leea30/git/fba.build/bubble/current/run_FlyDiscoRegisterTrx.sh /groups/branson/bransonlab/projects/olympiad/MCR/v717 /tier2/branson/fly_bubble/01_tracked/social_GMR_22D03_AE_01_dTrpA1_RigA_20160325T134540 
% settingsdir = ' /groups/flyprojects/home/olympiad/bubble_bin/FlyBubbleAnalysis/settings/ analysis_protocol current_bubble'
%%
modpath
%%
settingsdir = '/groups/rubin/home/schretterc/Documents/FlyDiscoAnalysis/FlyDiscoAnalysis-new/settings';
% rootdir = '/tier2/branson/FlyBubble_testing/social_flies';
expdir = '/groups/rubin/home/schretterc/Documents/FlyBowl_MingBowlTests/20XUAS_CsChrimson_mVenus_attP18_SS36564_flyBowlMing_20211211_continuous_Ramp10to90_20211212T120345';
%% parameters
% analysis_protocol = 'current_bubble';
% analysis_protocol = '20160415_flybubble_non_optogenetic';
% analysis_protocol = 'current_bubble';
% analysis_protcol = '20211014_flybubbleRed_LED';
analysis_protocol = '20211210_FlyBowlRED';
params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol};
%% run once
FlyDiscoDetectIndicatorLedOnOff(expdir,params{:});

