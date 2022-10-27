do_reset_destination_folder = true ;
do_use_bqueue = false ;
do_actually_submit_jobs = true ;

% Where does this script live?
this_script_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_script_path) ;
fly_disco_analysis_folder_path = fileparts(this_folder_path) ;
flydisco_folder_path = fileparts(fly_disco_analysis_folder_path) ;
root_example_experiments_folder_path = fullfile(flydisco_folder_path, 'example-experiments') ;
read_only_example_experiments_folder_path = fullfile(root_example_experiments_folder_path, 'single-passing-test-suite-experiment-read-only') ;

% Specify the "per-lab" configuration here
cluster_billing_account_name = 'scicompsoft' ;
rig_host_name = 'beet.hhmi.org' ;
rig_user_name = 'bransonk' ;
rig_data_folder_path = '/cygdrive/e/flydisco_data/scicompsoft' ;
goldblum_destination_folder_path = fullfile(root_example_experiments_folder_path, 'test-goldblum-destination-folder') ;
settings_folder_path = fullfile(fly_disco_analysis_folder_path, 'settings-internal') ;
per_lab_configuration = struct() ;
per_lab_configuration.cluster_billing_account_name = cluster_billing_account_name ;
per_lab_configuration.host_name_from_rig_index = {rig_host_name} ;
per_lab_configuration.rig_user_name_from_rig_index = {rig_user_name} ;
per_lab_configuration.data_folder_path_from_rig_index = {rig_data_folder_path} ;
per_lab_configuration.destination_folder = goldblum_destination_folder_path ;    
per_lab_configuration.settings_folder_path = settings_folder_path ;
%per_lab_configuration.does_use_per_user_folders = true ;

% Get the relative paths of all the experiment folders
absolute_path_to_read_only_folder_from_experiment_index = find_experiment_folders(read_only_example_experiments_folder_path) ;
relative_path_to_folder_from_experiment_index = ...
    cellfun(@(abs_path)(relpath(abs_path, read_only_example_experiments_folder_path)), ...
            absolute_path_to_read_only_folder_from_experiment_index, ...
            'UniformOutput', false) ;

if do_reset_destination_folder ,        
    % Delete the destination folder
    if exist(goldblum_destination_folder_path, 'file') ,
        return_code = system_from_list_with_error_handling({'rm', '-rf', goldblum_destination_folder_path}) ;
    end

    % % Recopy the analysis test folder from the template
    % fprintf('Resetting analysis test folder...\n') ;
    % read_only_example_experiments_folder_path = fullfile(root_example_experiments_folder_path, 'test-goldblum-example-experiments-folder') ;
    % reset_goldblum_example_experiments_working_copy_folder(read_only_example_experiments_folder_path, read_only_example_experiments_folder_path) ;

    % Copy to the destination folder
    rig_lab_data_folder_path = rig_data_folder_path ;
    fprintf('Transfering data to the destination path...\n') ;
    ensure_folder_exists(fileparts(goldblum_destination_folder_path)) ;  %#ok<UNRCH>
    command_line = {'cp', '-R', '-T', read_only_example_experiments_folder_path, goldblum_destination_folder_path} ;
    system_from_list_with_error_handling(command_line) ;
    % Should make goldblum_destination_folder_path a clone of
    % example_experiments_folder_path, since we use -T option
else
    % Just least goldblum_destination_folder_path as-is
end

% Add symlinks to the to-process folder so that they will actually get processed
folder_path_from_experiment_index = find_experiment_folders(goldblum_destination_folder_path) ;
to_process_folder_path = fullfile(goldblum_destination_folder_path, 'to-process') ;
ensure_folder_exists(to_process_folder_path) ;
experiment_count = length(folder_path_from_experiment_index) ;
for i = 1 : experiment_count ,
    experiment_folder_path = folder_path_from_experiment_index{i} ;
    [~, experiment_folder_name] = fileparts2(experiment_folder_path) ;
    symlink_path = fullfile(to_process_folder_path, experiment_folder_name) ;
    command_line = {'ln', '-s', '-f', '-T', experiment_folder_path, symlink_path} ;
    system_from_list_with_error_handling(command_line) ;
end

% Run goldblum
analysis_parameters = { 'doplotperframestats', 'on' } ;
fprintf('Running goldblum...\n') ;
do_transfer_data_from_rigs = false ;
do_run_analysis = true ;
goldblum(do_transfer_data_from_rigs, do_run_analysis, do_use_bqueue, do_actually_submit_jobs, analysis_parameters, per_lab_configuration) ;        

% Check that the expected files are present on dm11
local_verify(read_only_example_experiments_folder_path, goldblum_destination_folder_path) ;

% Check that some of the expected outputs were generated
all_tests_passed_from_experiment_index = check_for_pipeline_output_files(relative_path_to_folder_from_experiment_index, goldblum_destination_folder_path) ;
if all(all_tests_passed_from_experiment_index) ,
    fprintf('All experiment folder checks pass at 1st check, except those that were expected not to pass.\n') ;
else
    relative_path_to_folder_from_failed_experiment_index = ...
        relative_path_to_folder_from_experiment_index(~all_tests_passed_from_experiment_index)   %#ok<NOPTS,NASGU>   
    error('Some experiments had problems at 1st check') ;
end

% Run goldblum again, make sure nothing has changed
goldblum(do_transfer_data_from_rigs, do_run_analysis, do_use_bqueue, do_actually_submit_jobs, analysis_parameters, per_lab_configuration) ;        

% % Check that the expected files are present on dm11
% example_experiments_folder_destination_path = fullfile(destination_folder_path, 'taylora', 'analysis-test-folder') ;
% local_verify(example_experiments_folder_path, example_experiments_folder_destination_path) ;

% Check that some of the expected outputs were generated
all_tests_passed_from_experiment_index = check_for_pipeline_output_files(relative_path_to_folder_from_experiment_index, goldblum_destination_folder_path) ;
if all(all_tests_passed_from_experiment_index) ,
    fprintf('All experimental folder checks pass at 2nd check, except those that were expected not to pass.\n') ;
else
    relative_path_to_folder_from_failed_experiment_index = ...
        relative_path_to_folder_from_experiment_index(~all_tests_passed_from_experiment_index)  %#ok<NOPTS,NASGU>
    error('Some experiments had problems at 2nd check') ;
end

% If get here, all is well
[~, this_script_name] = fileparts(this_script_path) ;
fprintf('All tests in %s.m passed.\n', this_script_name) ;
