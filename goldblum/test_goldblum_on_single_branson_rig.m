% Test by copying some experiments to a single Branson Lab rig machine, then
% using goldblum to suck the data back and analyze the experiments

% Set some options
do_use_bqueue = true ;
do_actually_submit_jobs = true ;
do_run_analysis_in_debug_mode = true ;

% Figure out where this script lives in the filesystem
this_script_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_script_path) ;

% This is a folder with several "raw" experiment folders in it.
% We don't write to this folder, we only read from it.
analysis_test_template_folder_path = fullfile(this_folder_path, 'analysis-test-template') ;

% This stuff goes into the per-lab configuration that goldblum uses
cluster_billing_account_name = 'scicompsoft' ;
remote_host_name = 'beet.hhmi.org' ;
remote_host_name_from_rig_index = { remote_host_name } ;
remote_user_name = 'bransonk' ;
rig_user_name_from_rig_index = { remote_user_name } ;
remote_data_root_folder_path = '/cygdrive/e/flydisco_data' ;
data_folder_path_from_rig_index = {remote_data_root_folder_path} ;
destination_folder_path = fullfile(this_folder_path, 'goldblum-test-destination-folder') ;    
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

% Delete the destination folder so we're starting fresh
if exist(destination_folder_path, 'file') ,
    return_code = system_with_error_handling(sprintf('rm -rf %s', destination_folder_path)) ;
end

% Delete the remote lab data folder
remote_experiments_folder_path = remote_data_root_folder_path ;
escaped_remote_experiments_folder_path = escape_path_for_bash(remote_experiments_folder_path) ;
command_line = sprintf('ssh %s@%s rm -rf %s', remote_user_name, remote_host_name, escaped_remote_experiments_folder_path) ;
system_with_error_handling(command_line) ;

% Create the remote lab data folder
command_line = sprintf('ssh %s@%s mkdir %s', remote_user_name, remote_host_name, escaped_remote_experiments_folder_path) ;
system_with_error_handling(command_line) ;

% % Recopy the analysis test folder from the template
% reset_analysis_test_folder() ;

% Copy experiments to the rig computer
experiment_names = simple_dir(analysis_test_template_folder_path) ;
for i = 1 : length(experiment_names) ,
    experiment_name = experiment_names{i} ;
    local_experiment_folder_path = fullfile(analysis_test_template_folder_path, experiment_name) ;
    command_line = sprintf('scp -r %s %s@%s:%s', local_experiment_folder_path, remote_user_name, remote_host_name, remote_experiments_folder_path) ;
    system_with_error_handling(command_line) ;
end

% Run goldblum
goldblum(do_use_bqueue, do_actually_submit_jobs, do_run_analysis_in_debug_mode, per_lab_configuration) ;        

% Check that the expected files are present on dm11
local_verify(analysis_test_template_folder_path, destination_folder_path) ;

% Check that some of the expected outputs were generated
test_file_names = {'perframe' 'scores_AttemptedCopulation.mat' 'scoresBackup.mat' 'registered_trx.mat' 'wingtracking_results.mat'} ;
for i = 1 : length(test_file_names) ,
    test_file_name = test_file_names{i} ;
    test_file_path = ...
        fullfile(destination_folder_path, ...
                 'SS36564_20XUAS_CsChrimson_mVenus_attP18_flyBowlMing_20200227_Continuous_2min_5int_20200107_20200229T132141', ...
                 test_file_name) ;
    if ~exist(test_file_path, 'file') ,
        error('No output file at %s', test_file_path) ;
    end
end    

% Check that the rig experiments folder is empty now
entry_names = ...
    list_remote_dir(remote_user_name, remote_host_name, remote_experiments_folder_path) ;
if ~isempty(entry_names) ,
    error('Remote user folder %s:%s is not empty', remote_host_name, remote_experiments_folder_path) ;
end

% Run goldblum again, make sure nothing has changed
goldblum(do_use_bqueue, do_actually_submit_jobs, do_run_analysis_in_debug_mode, per_lab_configuration) ;        

% Check that the expected files are present on dm11
local_verify(analysis_test_template_folder_path, destination_folder_path) ;

% Check that some of the expected outputs were generated
test_file_names = {'perframe' 'scores_AttemptedCopulation.mat' 'scoresBackup.mat' 'registered_trx.mat' 'wingtracking_results.mat'} ;
for i = 1 : length(test_file_names) ,
    test_file_name = test_file_names{i} ;
    test_file_path = ...
        fullfile(destination_folder_path, ...
                 'SS36564_20XUAS_CsChrimson_mVenus_attP18_flyBowlMing_20200227_Continuous_2min_5int_20200107_20200229T132141', ...
                 test_file_name) ;
    if ~exist(test_file_path, 'file') ,
        error('No output file at %s', test_file_path) ;
    end
end    

% Check that the rig experiments folder is empty now
entry_names = ...
    list_remote_dir(remote_user_name, remote_host_name, remote_experiments_folder_path) ;
if ~isempty(entry_names) ,
    error('Remote user folder %s:%s is not empty', remote_host_name, remote_experiments_folder_path) ;
end

% If get here, all is well
[~, this_script_name] = fileparts(this_script_path) ;
fprintf('All tests in %s.m passed.\n', this_script_name) ;
