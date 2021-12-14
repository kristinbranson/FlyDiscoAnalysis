%DebugFlyDiscoDectectIndicatorLedOnOff.m 4/25/2016
% 
% /groups/flyprojects/home/leea30/git/fba.build/bubble/current/run_FlyDiscoRegisterTrx.sh /groups/branson/bransonlab/projects/olympiad/MCR/v717 /tier2/branson/fly_bubble/01_tracked/social_GMR_22D03_AE_01_dTrpA1_RigA_20160325T134540 
% settingsdir = ' /groups/flyprojects/home/olympiad/bubble_bin/FlyBubbleAnalysis/settings/ analysis_protocol current_bubble'
% %%
% addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
% addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
% settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings';
% %rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
modpath
%%
settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/settings';


% rootdir = '/tier2/branson/FlyBubble_testing/social_flies';
% expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20211014_testing_FlyBubbleRed/socialCsChr_GMR_72C11_AE_01_CsChrimson_RigD_20191114T172654';
% expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20211129_testingLEDdetectioncount/NorpA_JHS_K_85321_RigD_20210902T075331';
expdirs = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20211129_testingLEDdetectioncount/VNC_JRC_SS50051_RigD_20210426T151519',...
    '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20211129_testingLEDdetectioncount/NorpA_JHS_K_85321_RigD_20210902T075331'};

%% parameters

% analysis_protocol = 'current_bubble';
% analysis_protocol = '20160415_flybubble_non_optogenetic';
% analysis_protocol = 'current_bubble';
% analysis_protcol = '20211014_flybubbleRed_LED';
analysis_protocol = '20210531_flybubble_LED_AR_20210819'; %has fix for Kristin's processing of Katie's data
params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol};

%% run once
for i = 1:numel(expdirs),
    expdir = expdirs{i};
    disp(expdir);
    try
        FlyDiscoDectectIndicatorLedOnOff(expdir,params{:});
    catch me
    end
end

