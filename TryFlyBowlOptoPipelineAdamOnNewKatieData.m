%Debug_FlyBubblePipeline.m

path_to_this_script = mfilename('fullpath') ;
% Set up path to libraries
modpath() ;

% Get the path to the settings folder
%settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings';
%settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis_github/settings';
path_to_this_folder = fileparts( mfilename('fullpath') ) ;
path_to_parent_folder = fileparts(path_to_this_folder) ;
settingsdir = fullfile(path_to_this_folder, 'settings') 

%% parameters
% analysis_protocol = '20161215_flybubble_seedsa';
% analysis_protocol = 'current_non_olympiad_rubin_grooming';
% analysis_protocol = 'current_bubble';
% analysis_protocol = '20190826_flybubble_pipelinemore';

%analysis_protocol = '20190712_flybubble_flybowloptoKatie_mingrig_avi';
analysis_protocol = '20190712_flybubble_flybowloptoKatie_mingrig'

% analysis_protocol = '20200123_flybubble_centralcomplex';
% analysis_protocol = '20150717_flybubble_flybowlMing';

params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol, ...
  'forcecompute',false,...
  'doautomaticchecksincoming',false,...
  'doregistration',true,...
  'doledonoffdetection',true,...  
  'dotrackwings',true,...
  'dosexclassification',true,...
  'docomputeperframefeatures',true,...
  'docomputehoghofperframefeatures',false,...
  'dojaabadetect',true,...
  'docomputeperframestats',false,...
  'doplotperframestats',false,...
  'domakectraxresultsmovie',true,...
  'doextradiagnostics',false,...
  'doanalysisprotocol',isunix,...
  'doautomaticcheckscomplete',false};

%% experiment directory
% expdir = '/groups/branson/home/robiea/Projects_data/temp/testing_pipelinescript';
% expdir = '/groups/branson/home/robiea/Projects_data/FlyBowl_Opto/TestData/bubble/cx_JHS_K_85321_CsChr_RigB_20150707T134343';
% expdir = '/groups/branson/home/robiea/Projects_data/Katie/tracked/20190716T130615_rig1_flyBowl3__GMR_OL0077B_rubin_protocol_OL0077_testing_noPixels';
% expdir = '/groups/branson/home/robiea/Projects_data/FlyBubble/bubble_data_social/testing_pipeline/social_GMR_71G01_AE_01_dTrpA1_RigD_20190808T165554';
% expdir = '/groups/branson/home/robiea/Projects_data/FlyBubble/bubble_data_social/testing_pipeline/social_GMR_70G01_AE_01_dTrpA1_RigC_20190808T165559';
% expdir = '/groups/branson/home/robiea/Projects_data/FlyBowl_Opto/TestData/NewOptoBowls/20191201T180528_rig1_flyBowl3__SS47478_20XUASCsChrimsonattp18_protocol_OL0077_testing_long';
% expdir = '/groups/branson/home/robiea/Projects_data/FlyBowl_Opto/TestData/NewOptoBowls/20191204T151156_rig1_flyBowl3__SS47478_20XUASCsChrimsonattp18_protocol_OL0077_testing_long_fortesting';
% expdir = '/groups/branson/home/robiea/Projects_data/FlyBowl_Opto/TestData/NewOptoBowls/20191204T151156_rig1_flyBowl3__SS47478_20XUASCsChrimsonattp18_protocol_OL0077_testing_long_fortesting_20191211';
% expdir = '/groups/branson/home/robiea/Projects_data/FlyBowl_Opto/TestData/NewOptoBowls/20191204T151156_rig1_flyBowl3__SS47478_20XUASCsChrimsonattp18_protocol_OL0077_testing_long_fortesting_TESTwingtracking';
% expdir = '/groups/branson/home/robiea/Projects_data/FlyBowl_Opto/TestData/NewOptoBowls/20191204T151156_rig1_flyBowl3__SS47478_20XUASCsChrimsonattp18_protocol_OL0077_testing_long_fortesting_20191211_TESTINGwingtracking';
% not finding indicator - timing is off
% expdir = '/groups/branson/home/robiea/Projects_data/FlyBowl_Opto/TestData/mingbowl/testing/emptysplit_20xUAS-ChrimsonRmVenusattp18_flyBowlMing_nopause_lengthofpersis_2min_10int_20191218T093239';
%expdir = '/groups/branson/home/robiea/Projects_data/FlyBowl_Opto/TestData/mingbowl/testing/emptysplit_20xUAS-ChrimsonRmVenusattp18_flyBowlMing_nopause_lengthofpersis_2min_10int_20191218T093239_2' ;

expdir = ...
    fullfile(path_to_parent_folder, ...
             'experiments', ...
             '20201014T102140_rig1_flyBowl1__2135003_1xLwt_attp40_4stop1_protocolRGBTemplate') 

%expdir =
%'/groups/branson/home/robiea/Projects_data/Katie/VisualOrientation/reprocessedbyAlice/SS32237_20xUASCsChrimsonDARK_flyBowlMing_20200227_Continuous_2min_5int_20200107_20200404T080309'
% ^ doesn't seem to exist

% calc hof/hog features
% expdir = '/groups/branson/home/robiea/Projects_data/Labeler_APT/socialCsChr_JRC_SS36551_CsChrimson_RigB_20190910T160839';
% expdir = '/groups/branson/bransonlab/DataforAPT/FlyBubble/Grooming/cx_26B12x26A02_CsChr_None_RigC_20150731T135550';
% expdir = '/groups/branson/bransonlab/DataforAPT/FlyBubble/Grooming/cx_26B12x24C08_CsChr_None_RigA_20150731T134907';
% expdir = '/groups/branson/bransonlab/DataforAPT/FlyBubble/Grooming/cx_60E02x52F12_CsChr_None_RigD_20150731T135953';
% expdir = '/groups/branson/bransonlab/DataforAPT/FlyBubble/Grooming/cx_76F12x18C04_CsChr_None_RigB_20150731T135428';
% testing hog/hof computation on linux
% expdir = '/groups/branson/bransonlab/DataforAPT/FlyBubble/Grooming/HOGHOF_inGUI/cx_26B12x26A02_CsChr_None_RigC_20150731T135550';

% explist = {'/groups/branson/home/robiea/Projects_data/Labeler_APT/socialCsChr_JRC_SS36551_CsChrimson_RigA_20190910T160809', ...
%     '/groups/branson/home/robiea/Projects_data/Labeler_APT/socialCsChr_JRC_SS56987_CsChrimson_RigB_20190910T163328'};
% expdir = '/groups/branson/home/robiea/Projects_data/Labeler_APT/socialCsChr_JRC_SS56987_CsChrimson_RigA_20190910T163337';
% not finding LED indicator - pulse IR LED indicator (delayed)
% expdir = '/groups/branson/home/robiea/Projects_data/FlyBowl_Opto/TestData/mingbowl/testing/SS36564_20XUAS_CsChrimson_mVenus_attP18_flyBowlMing_Pulse_5min_5int_20200203_20200205T072836';
% expdir = '/groups/branson/home/robiea/Projects_data/FlyBowl_Opto/TestData/mingbowl/testing/SS56987_20XUAS_CsChrimson_mVenus_attP18_flyBowlMing_Pulse_5min_5int_20200203_20200205T071450';
% fixed by becca
% expdir = '/groups/branson/home/robiea/Projects_data/FlyBowl_Opto/TestData/mingbowl/testing/emptysplit_20XUAS_CsChrimson_mVenus_attP18_flyBowlMing_Continuous_2min_10int_20200107_20200114T092950';
%% run
% for i = 1:numel(explist)
%     expdir = explist{i};
% [success,msgs,stage] = FlyBubblePipeline(expdir,params{:});
% end

[success,msgs,stage] = FlyBubblePipeline(expdir,params{:});
