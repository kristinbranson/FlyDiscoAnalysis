%% set path
modpath

                                                    
%% set stuff up
cluster_billing_account_name = 'rubin' ;
% only used for setting configuration of FlyDiscoPipeline 
user_name_for_configuration_purposes = 'rubinlab' ;  
do_use_bqueue = false ;
do_actually_submit_jobs = false ;
ssh_host_name = 'submit.int.janelia.org' ;


settings_folder_path = '/groups/rubin/home/rubinlabuser/RubinFlyDiscoSettings/RubinFlyDiscoSettings/settings';
analysis_protocol = '20230302_BubblefRGB_GtAg';

optional_argument_list = ...
    {'settingsdir', settings_folder_path, ...
    'analysis_protocol',analysis_protocol, ...
     'do_try', false} ; 

%% experiments 

folder_path_from_experiment_index = {'/groups/rubin/home/schretterc/Documents/FlyBubblefRGB_Chrim/ForTestingAPT/20230222T161958_rig1_flyBubble1_SS37072_20XUASCsChrimsonattp18_20221129_AggAVLPChrim'};

%% Call the testing function to do the real work
test_transfero_FlyDiscoPipeline_wrapper_on_experiment_list(folder_path_from_experiment_index, ...
                                                           cluster_billing_account_name, ...
                                                           user_name_for_configuration_purposes, ...
                                                           do_use_bqueue, ...
                                                           do_actually_submit_jobs, ...
                                                           ssh_host_name, ...
                                                           optional_argument_list{:}) 
