
%% set path
modpath

                                                    
%% set stuff up
cluster_billing_account_name = 'branson' ;
% only used for setting configuration of FlyDiscoPipeline 
user_name_for_configuration_purposes = 'bransonlab' ;  
do_use_bqueue = false ;
do_actually_submit_jobs = false ;
ssh_host_name = 'login1.int.janelia.org' ;


settings_folder_path = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/settings-internal';
% analysis_protocol = '20230206_BubblefRGB_GtAg';
% VNC2 experiments
analysis_protocol = '20230213_flybubble_LED_VNC2';
% VNC experiments
%analysis_protocol = '20230308_flybubble_LED_VNC';

% optional_argument_list = ...
%     {'settingsdir', settings_folder_path, ...
%     'analysis_protocol',analysis_protocol, ...
%      'do_try', false} ; 
optional_argument_list = ...
    {'settingsdir', settings_folder_path, ... 
     'do_try', false} ;

% optional_argument_list = {'do_try', false} ;

%% experiments 

% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20230206_testAPTtrackingPipeline/20230131T080157_rig2_flyBubble1_SSemptyCtoDex_20XUASCsChrimsonattp18_20221129_AggAVLPChrim'};
% folder_path_from_experiment_index ={'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20230206_testAPTtrackingPipeline/20230202T084401_rig2_flyBubble1_SS37072Dex_20XUASCsChrimsonattp18_20230110_AggAVLPChrim96882'};
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20230201_movie4Thursdayseminar/VNC2_YNA_K_162984_RigA_20220518T115615'};
% testing entire pipeline
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20230206_testAPTtrackingPipeline/VNC2_YNA_K_162984_RigA_20220518T115615'};
% testing perframe, JAABA, APT, aptresultsmovie
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20230215_testingnewStage/VNC2_JRC_SS71988_RigC_20220921T115817'};
% testing APT merge (checked out current APT commit on local submodule branch)
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20230206_testAPTtrackingPipeline/VNC2_JRC_SS71988_RigB_20220921T115713'}
% trying to reprocess Roman's training data movie 
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20230306_RomanTrainingDataSet/VNC_YNA_K_162984_RigD_20210526T155055_fromdiscodata'};
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20230509_autochecksERROR/VNC2_EXT_VGLUT-GAL4_RigD_20230502T112921'};
folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20230509_autochecksERROR/VNC2_YNA_K_162984_RigD_20230503T135414'};
%% Call the testing function to do the real work
test_transfero_FlyDiscoPipeline_wrapper_on_experiment_list(folder_path_from_experiment_index, ...
                                                           cluster_billing_account_name, ...
                                                           user_name_for_configuration_purposes, ...
                                                           do_use_bqueue, ...
                                                           do_actually_submit_jobs, ...
                                                           ssh_host_name, ...
                                                           optional_argument_list{:}) 