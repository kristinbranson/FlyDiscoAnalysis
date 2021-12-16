% Test by copying some experiments to a single Branson Lab rig machine, then
% using goldblum to suck the data back and analyze the experiments

% Set some options
do_use_bqueue = true ;
do_actually_submit_jobs = true ;
do_run_analysis_in_debug_mode = false ;

% Figure out where this script lives in the filesystem
this_script_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_script_path) ;

% This is a folder with several "raw" experiment folders in it.
% We don't write to this folder, we only read from it.
analysis_test_template_folder_path = fullfile(this_folder_path, 'analysis-test-template') ;

% This stuff goes into the per-lab configuration that goldblum uses
cluster_billing_account_name = 'branson' ;
remote_host_name = 'beet.hhmi.org' ;
remote_host_name_from_rig_index = { remote_host_name } ;
remote_user_name = 'bransonk' ;
rig_user_name_from_rig_index = { remote_user_name } ;
remote_data_root_folder_path = '/cygdrive/e/flydisco_data' ;
data_folder_path_from_rig_index = {remote_data_root_folder_path} ;
destination_folder_path = '/groups/branson/bransonlab/flydisco_data' ;
settings_folder_path = fullfile(this_folder_path, 'FlyDiscoAnalysis/settings') ;
does_use_per_user_folders = false ;

% Specify the "per-lab" configuration here
per_lab_configuration = struct() ;
per_lab_configuration.cluster_billing_account_name = cluster_billing_account_name ;
per_lab_configuration.host_name_from_rig_index = remote_host_name_from_rig_index ;
per_lab_configuration.rig_user_name_from_rig_index = rig_user_name_from_rig_index ;
per_lab_configuration.data_folder_path_from_rig_index = data_folder_path_from_rig_index ;
per_lab_configuration.destination_folder = destination_folder_path ;    
per_lab_configuration.settings_folder_path = settings_folder_path ;
per_lab_configuration.does_use_per_user_folders = does_use_per_user_folders ;

% Run goldblum
goldblum(do_use_bqueue, do_actually_submit_jobs, do_run_analysis_in_debug_mode, per_lab_configuration) ;        

