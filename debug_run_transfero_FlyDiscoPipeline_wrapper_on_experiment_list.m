% run pipeline locally AR 20231116


cluster_billing_account_name = 'branson' ;
user_name_for_configuration_purposes = 'bransonlab' ; 

settings_folder_path = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings';
optional_argument_list = {'settingsdir', settings_folder_path};
do_use_bqueue = false ;
do_actually_submit_jobs = false ;
do_try = false ;
ssh_host_name = 'login2.int.janelia.org' ;

% explist to run
folder_path_from_experiment_index = {
    '/groups/branson/bransonlab/alice/LPC1test/LPC1ShiKir_BJD_SS02407_RigA_20231116T150618'};
% '/groups/branson/bransonlab/alice/LPC1test/LPC1CsChr_BJD_SS02407_RigA_20231116T132456'}; % started OK


% Call the user-facing run function to do the real work
run_transfero_FlyDiscoPipeline_wrapper_on_experiment_list(folder_path_from_experiment_index, ...
                                                          cluster_billing_account_name, ...
                                                          user_name_for_configuration_purposes, ...
                                                          do_use_bqueue, ...
                                                          do_actually_submit_jobs, ...
                                                          do_try, ...
                                                          ssh_host_name, ...
                                                          optional_argument_list{:}) 