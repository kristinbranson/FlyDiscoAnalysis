% run pipeline locally AR 20231116


cluster_billing_account_name = 'branson' ;
user_name_for_configuration_purposes = 'bransonlab' ; 

% settings_folder_path = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings';
settings_folder_path = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/settings-internal';
analysis_protocol_use = '20240124_multibubble_secondtry';
optional_argument_list = {'settingsdir', settings_folder_path,'analysis_protocol',analysis_protocol_use};
do_use_bqueue = false ;
do_actually_submit_jobs = false ;
do_try = false ;
ssh_host_name = 'login2.int.janelia.org' ;

% explist to run
% folder_path_from_experiment_index = {
%     '/groups/branson/bransonlab/alice/LPC1test/LPC1ShiKir_BJD_SS02407_RigA_20231116T150618'};
% '/groups/branson/bransonlab/alice/LPC1test/LPC1CsChr_BJD_SS02407_RigA_20231116T132456'}; % started OK

% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20231117_testCleanOut/LPC1ShiKir_BJD_SS02407_RigA_20231116T150618'};
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20231218_testingsettingsdirsDNs/DNcontrol_GMR_SS02575_RigB_20231213T135549'}
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20231218_testingsettingsdirsDNs/DNretinal_SS77821_RigA_20231213T141155'}
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/multibubble_data/20240124_testingpipeline/20240122T154843_rig2_flyBowl1_1xLwt_attp40_4stop1_3xDSCPLwt_attp40_3stop1_SlowRamp'}
folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/multibubble_data/20240124_testingpipeline/NewMultiBubble_14Feb2024/20240214T113854_rig1_multibubble_71G01_3_1_noCS_SlowRamp'};
% '/groups/branson/home/robiea/Projects_data/FlyDisco/multibubble_data/20240124_testingpipeline/NewMultiBubble_14Feb2024/20240214T112101_rig1_multibubble_71G01_1_1_noCS_SlowRamp'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/multibubble_data/20240124_testingpipeline/NewMultiBubble_14Feb2024/20240214T110436_rig1_multibubble_71G01_2_1_noCS_SlowRamp'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/multibubble_data/20240124_testingpipeline/NewMultiBubble_14Feb2024/20240214T104810_rig1_multibubble_71G01_1_2_noCS_SlowRamp'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/multibubble_data/20240124_testingpipeline/NewMultiBubble_14Feb2024/20240214T102958_rig1_multibubble_CSx71G01_1_2_UAS_Chrimson_Venus_X_0070_SlowRamp'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/multibubble_data/20240124_testingpipeline/NewMultiBubble_14Feb2024/20240214T101329_rig1_multibubble_CSx71G01_2_1_UAS_Chrimson_Venus_X_0070_SlowRamp'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/multibubble_data/20240124_testingpipeline/NewMultiBubble_14Feb2024/20240214T095625_rig1_multibubble_CSx71G01_1_1_UAS_Chrimson_Venus_X_0070_SlowRamp'
% '/groups/branson/home/robiea/Projects_data/FlyDisco/multibubble_data/20240124_testingpipeline/NewMultiBubble_14Feb2024/20240214T093136_rig1_multibubble_CSx71G01_3_1_UAS_Chrimson_Venus_X_0070_SlowRamp'}
% Call the user-facing run function to do the real work
run_transfero_FlyDiscoPipeline_wrapper_on_experiment_list(folder_path_from_experiment_index, ...
                                                          cluster_billing_account_name, ...
                                                          user_name_for_configuration_purposes, ...
                                                          do_use_bqueue, ...
                                                          do_actually_submit_jobs, ...
                                                          do_try, ...
                                                          ssh_host_name, ...
                                                          optional_argument_list{:}) 