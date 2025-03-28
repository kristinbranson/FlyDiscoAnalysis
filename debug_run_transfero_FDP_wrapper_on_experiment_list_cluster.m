% % run pipeline locally AR 20231116


cluster_billing_account_name = 'branson' ;
user_name_for_configuration_purposes = 'bransonlab' ; 

% settings_folder_path = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings';
% settings_folder_path = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/settings-internal';
% settings_folder_path = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings';
% analysis_protocol_use = '20240124_multibubble_secondtry';
% analysis_protocol_use = '20240521_flybubble_LED_VNC3';
% analysis_protocol_use   = '20240908_flybubble_centralcomplex';
% analysis_protocol_use   = '20240912_flybubble_socialCsChr';
% analysis_protocol_use = '20240910_flybubblered_nochr_flytracker';
% analysis_protocol_use = '20240917_flybubblered_nochr_flytracker_finaljab';
% analysis_protocol_use = '20241203_flybubble_LED_VNC3';
settings_folder_path = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings';
analysis_protocol_use =  '20250111_flybubble_LED_NorpA';

optional_argument_list = {'settingsdir', settings_folder_path,'analysis_protocol',analysis_protocol_use};
do_use_bqueue = false ;
do_actually_submit_jobs = false ;
do_try = true ;
ssh_host_name = 'login2.int.janelia.org' ;

% explist to run
% folder_path_from_experiment_index = {
%     '/groups/branson/bransonlab/alice/LPC1test/LPC1ShiKir_BJD_SS02407_RigA_20231116T150618'};
% '/groups/branson/bransonlab/alice/LPC1test/LPC1CsChr_BJD_SS02407_RigA_20231116T132456'}; % started OK

% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20231117_testCleanOut/LPC1ShiKir_BJD_SS02407_RigA_20231116T150618'};
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20231218_testingsettingsdirsDNs/DNcontrol_GMR_SS02575_RigB_20231213T135549'}
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20231218_testingsettingsdirsDNs/DNretinal_SS77821_RigA_20231213T141155'}
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/multibubble_data/20240124_testingpipeline/20240122T154843_rig2_flyBowl1_1xLwt_attp40_4stop1_3xDSCPLwt_attp40_3stop1_SlowRamp'}
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/multibubble_data/20240124_testingpipeline/NewMultiBubble_14Feb2024/20240214T113854_rig1_multibubble_71G01_3_1_noCS_SlowRamp'};
% '/groups/branson/home/robiea/Projects_data/FlyDisco/multibubble_data/20240124_testingpipeline/NewMultiBubble_14Feb2024/20240214T112101_rig1_multibubble_71G01_1_1_noCS_SlowRamp'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/multibubble_data/20240124_testingpipeline/NewMultiBubble_14Feb2024/20240214T110436_rig1_multibubble_71G01_2_1_noCS_SlowRamp'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/multibubble_data/20240124_testingpipeline/NewMultiBubble_14Feb2024/20240214T104810_rig1_multibubble_71G01_1_2_noCS_SlowRamp'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/multibubble_data/20240124_testingpipeline/NewMultiBubble_14Feb2024/20240214T102958_rig1_multibubble_CSx71G01_1_2_UAS_Chrimson_Venus_X_0070_SlowRamp'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/multibubble_data/20240124_testingpipeline/NewMultiBubble_14Feb2024/20240214T101329_rig1_multibubble_CSx71G01_2_1_UAS_Chrimson_Venus_X_0070_SlowRamp'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/multibubble_data/20240124_testingpipeline/NewMultiBubble_14Feb2024/20240214T095625_rig1_multibubble_CSx71G01_1_1_UAS_Chrimson_Venus_X_0070_SlowRamp'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/multibubble_data/20240124_testingpipeline/NewMultiBubble_14Feb2024/20240214T093136_rig1_multibubble_CSx71G01_3_1_UAS_Chrimson_Venus_X_0070_SlowRamp'}

% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240508_VNC3/VNC3_YNA_K_162984_RigA_20240507T182023'};
%test apt features in compute perframe and apt pff computation
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240521_testAPTfeatures/VNC2_YNA_K_162984_RigB_20220831T124607'};

% testing new stage
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset_inpipeline/VNC_YNA_K_162984_RigD_20210922T133035'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset_inpipeline/VNC_YNA_K_162984_RigC_20210922T132938'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset_inpipeline/VNC_YNA_K_162984_RigB_20210922T132823'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset_inpipeline/VNC_YNA_K_162984_RigA_20210922T132731'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset_inpipeline/VNC_YNA_K_162984_RigD_20210921T140825'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset_inpipeline/VNC_YNA_K_162984_RigC_20210921T140740'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset_inpipeline/VNC_YNA_K_162984_RigB_20210921T140516'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset_inpipeline/VNC_YNA_K_162984_RigA_20210921T140433'};

% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC2_YNA_K_162984_RigA_20220802T121645'}
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC2_YNA_K_162984_RigA_20220803T112613'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC2_YNA_K_162984_RigB_20220803T112705'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC2_YNA_K_162984_RigC_20220802T121835'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC2_YNA_K_162984_RigC_20220803T112800'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC2_YNA_K_162984_RigD_20220802T121919'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240611_testingboutmetrics_ctrlset/VNC2_YNA_K_162984_RigD_20220803T112857'};
% folder_path_from_experiment_index = {'/groups/branson/bransonlab/flydisco_data/VNC2_YNA_K_162984_RigA_20220802T121645'};

%test rerunning failed experiments
% folder_path_from_experiment_index = {'/groups/branson/bransonlab/flydisco_data/VNC3_YNA_K_162984_RigD_20240702T123257'};

%test new settings dir for cx
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyBubbleMethods/Data/CsChr_OL0046B/cx_GMR_OL0046B_CsChr_RigC_20150408T142551'};
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyBubbleMethods/Data/CsChr_OL0046B/cx_GMR_OL0046B_CsChr_RigA_20150617T094711'
% '/groups/branson/home/robiea/Projects_data/FlyBubbleMethods/Data/CsChr_OL0046B/cx_GMR_OL0046B_CsChr_RigB_20150623T102758'
% '/groups/branson/home/robiea/Projects_data/FlyBubbleMethods/Data/CsChr_OL0046B/cx_GMR_OL0046B_CsChr_RigD_20150708T092748'}

% test new settings for nochr
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyBubbleMethods/Data/JAABA/nochr_TrpA71G01_Unknown_RigC_20201216T163316'};
%run other 2 for JAABA
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyBubbleMethods/Data/JAABA/nochr_TrpA91B01_Unknown_RigB_20201212T170551'
% '/groups/branson/home/robiea/Projects_data/FlyBubbleMethods/Data/JAABA/nochr_TrpApBDP_Unknown_RigB_20201212T165050'};

% %46B expdirs
% folder_path_from_experiment_index = {'//groups/branson/home/robiea/Projects_data/FlyBubbleMethods/Data/CsChr_OL0046B/cx_GMR_OL0046B_CsChr_RigA_20150707T132832'
% '/groups/branson/home/robiea/Projects_data/FlyBubbleMethods/Data/CsChr_OL0046B/cx_GMR_OL0046B_CsChr_RigB_20150708T092121'
% '/groups/branson/home/robiea/Projects_data/FlyBubbleMethods/Data/CsChr_OL0046B/cx_GMR_OL0046B_CsChr_RigC_20150707T093155'
% '/groups/branson/home/robiea/Projects_data/FlyBubbleMethods/Data/CsChr_OL0046B/cx_GMR_OL0046B_CsChr_RigD_20150708T092748'};

% socialcsChr with cx settings dir 
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyBubbleMethods/Data/JAABA/socialCsChr_JRC_SS36551_CsChrimson_RigB_20190910T160839'};

% test experiemnts for social touch
% folder_path_from_experiment_index = textread('/groups/branson/home/robiea/Projects_data/FlyBubbleMethods/Data/TrpA/explist.txt','%s');
%replacement experiments for failed in list above
% folder_path_from_experiment_index =  {'/groups/branson/home/robiea/Projects_data/FlyBubbleMethods/Data/TrpA/nochr_TrpA71G01_Unknown_RigA_20201216T153505a'
% '/groups/branson/home/robiea/Projects_data/FlyBubbleMethods/Data/TrpA/nochr_TrpA91B01_Unknown_RigA_20201216T152414a'};

% social touch test list
% folder_path_from_experiment_index = textread('/groups/branson/home/robiea/Projects_data/FlyBubbleMethods/Data/TrpA/explist.txt','%s');
% test locomotion stage
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20241203_testinglocomotion/VNC3_JRC_SS100086_RigD_20240820T125028'};
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/Labeler_APT/FlyTracker/cx_JHS_K_85321_CsChr_RigD_20150909T163219'
% '/groups/branson/home/robiea/Projects_data/Labeler_APT/FlyTracker/cx_GMR_SS00168_CsChr_RigD_20150909T111218'
% '/groups/branson/home/robiea/Projects_data/Labeler_APT/FlyTracker/cx_GMR_SS00038_CsChr_RigB_20150729T150617'};

% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/Labeler_APT/FlyTracker/cx_GMR_SS00038_CsChr_RigB_20150729T150617'};

% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/Labeler_APT/FlyTracker/to_reprocess/nochr_TrpA20A02_Unknown_RigB_20201216T154643'
% '/groups/branson/home/robiea/Projects_data/Labeler_APT/FlyTracker/to_reprocess/nochr_TrpA71G01_Unknown_RigD_20201216T153831'
% '/groups/branson/home/robiea/Projects_data/Labeler_APT/FlyTracker/to_reprocess/socialCsChr_JRC_SS56987_CsChrimson_RigB_20190910T163328'};

folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20250111_testingFlybubble_social/NorpA5_JRC_SS56987_RigC_20210923T085517'};

% Call the user-facing run function to do the real work
run_transfero_FlyDiscoPipeline_wrapper_on_experiment_list(folder_path_from_experiment_index, ...
                                                          cluster_billing_account_name, ...
                                                          user_name_for_configuration_purposes, ...
                                                          do_use_bqueue, ...
                                                          do_actually_submit_jobs, ...
                                                          do_try, ...
                                                          ssh_host_name, ...
                                                          optional_argument_list{:}) 