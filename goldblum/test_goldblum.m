do_transfer_data_from_rigs = false ;
do_run_analysis = true ;
do_use_bqueue = true ;
do_actually_submit_jobs = true ;

% Where does this script live?
this_script_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_script_path) ;
fly_disco_analysis_folder_path = fileparts(this_folder_path) ;

% Specify the "per-lab" configuration here
lab_head_last_name = 'scicompsoft' ;
rig_host_name = 'arrowroot.hhmi.org' ;
rig_user_name = 'bransonk' ;
rig_data_folder_path = '/cygdrive/e/flydisco_data' ;
goldblum_destination_folder_path = fullfile(this_folder_path, 'goldblum-test-destination-folder') ;
settings_folder_path = fullfile(fly_disco_analysis_folder_path, 'settings') ;
per_lab_configuration = struct() ;
per_lab_configuration.lab_head_last_name = lab_head_last_name ;
per_lab_configuration.host_name_from_rig_index = {rig_host_name} ;
per_lab_configuration.rig_user_name_from_rig_index = {rig_user_name} ;
per_lab_configuration.data_folder_path_from_rig_index = {rig_data_folder_path} ;
per_lab_configuration.destination_folder = goldblum_destination_folder_path ;    
per_lab_configuration.settings_folder_path = settings_folder_path ;
%per_lab_configuration.does_use_per_user_folders = true ;

% Delete the destination folder
if exist(goldblum_destination_folder_path, 'file') ,
    return_code = system_with_error_handling(sprintf('rm -rf %s', goldblum_destination_folder_path)) ;
end

% Recopy the analysis test folder from the template
fprintf('Resetting analysis test folder...\n') ;
example_experiments_folder_path = reset_goldblum_example_experiments_working_copy_folder() ;

% Copy it to the rig computer (or just direct to the destination folder)
rig_lab_data_folder_path = fullfile(rig_data_folder_path, lab_head_last_name) ;
if do_transfer_data_from_rigs ,
    fprintf('Transfering data to the rig computer...\n') ;  %#ok<UNRCH>
    command_line = sprintf('scp -B -r %s/* %s@%s:%s', example_experiments_folder_path, rig_user_name, rig_host_name, rig_lab_data_folder_path) ; %#ok<UNRCH>
    system_with_error_handling(command_line) ;
else
    fprintf('Transfering data to the destination path...\n') ;
    ensure_folder_exists(goldblum_destination_folder_path) ;  %#ok<UNRCH>
    command_line = sprintf('cp -R %s/* %s', example_experiments_folder_path, goldblum_destination_folder_path) ;
    system_with_error_handling(command_line) ;
    
    % Add symlinks to the to-process folder so that they will actually get processed
    folder_path_from_experiment_index = find_experiment_folders(goldblum_destination_folder_path) ;
    to_process_folder_path = fullfile(goldblum_destination_folder_path, 'to-process') ;
    escaped_to_process_folder_path = escape_string_for_bash(to_process_folder_path) ;
    ensure_folder_exists(to_process_folder_path) ;
    experiment_count = length(folder_path_from_experiment_index) ;
    for i = 1 : experiment_count ,
        experiment_folder_path = folder_path_from_experiment_index{i} ;
        escaped_experiment_folder_path = escape_string_for_bash(experiment_folder_path) ;
        command_line = sprintf('ln -s %s %s', escaped_experiment_folder_path, escaped_to_process_folder_path) ;
        system_with_error_handling(command_line) ;        
    end    
end

% Run goldblum
fprintf('Running goldblum...\n') ;
goldblum(do_transfer_data_from_rigs, do_run_analysis, do_use_bqueue, do_actually_submit_jobs, [], per_lab_configuration) ;        

% Check that the expected files are present on dm11
local_verify(example_experiments_folder_path, goldblum_destination_folder_path) ;

% Check that some of the expected outputs were generated
test_file_names = {'perframe' 'scores_AttemptedCopulation.mat' 'scoresBackup.mat' 'registered_trx.mat' 'wingtracking_results.mat'} ;
for i = 1 : length(test_file_names) ,
    test_file_name = test_file_names{i} ;
    test_file_path = ...
        fullfile(goldblum_destination_folder_path, ...
                 'SS36564_20XUAS_CsChrimson_mVenus_attP18_flyBowlMing_20200227_Continuous_2min_5int_20200107_20200229T132141', ...
                 test_file_name) ;
    if ~exist(test_file_path, 'file') ,
        error('No output file at %s', test_file_path) ;
    end
end    

% Check that the rig lab folder is empty now
if do_transfer_data_from_rigs ,
    entry_names = ...
        list_remote_dir(rig_user_name, rig_host_name, rig_lab_data_folder_path) ;  %#ok<UNRCH>
    if ~isempty(entry_names) ,
        error('Rig lab data folder %s:%s is not empty', rig_host_name, rig_lab_data_folder_path) ;
    end
end

% Run goldblum again, make sure nothing has changed
goldblum(do_transfer_data_from_rigs, do_run_analysis, do_use_bqueue, do_actually_submit_jobs, [], per_lab_configuration) ;        

% % Check that the expected files are present on dm11
% example_experiments_folder_destination_path = fullfile(destination_folder_path, 'taylora', 'analysis-test-folder') ;
% local_verify(example_experiments_folder_path, example_experiments_folder_destination_path) ;

% Check that some of the expected outputs were generated
test_file_names = {'perframe' 'scores_AttemptedCopulation.mat' 'scoresBackup.mat' 'registered_trx.mat' 'wingtracking_results.mat'} ;
for i = 1 : length(test_file_names) ,
    test_file_name = test_file_names{i} ;
    test_file_path = ...
        fullfile(goldblum_destination_folder_path, ...
                 'SS36564_20XUAS_CsChrimson_mVenus_attP18_flyBowlMing_20200227_Continuous_2min_5int_20200107_20200229T132141', ...
                 test_file_name) ;
    if ~exist(test_file_path, 'file') ,
        error('No output file at %s', test_file_path) ;
    end
end    

% Check that the rig lab folder is empty now
if do_transfer_data_from_rigs ,
    entry_names = ...
        list_remote_dir(rig_user_name, rig_host_name, rig_lab_data_folder_path) ;  %#ok<UNRCH>
    if ~isempty(entry_names) ,
        error('Remote user folder %s:%s is not empty', rig_host_name, rig_lab_data_folder_path) ;
    end
end

% If get here, all is well
[~, this_script_name] = fileparts(this_script_path) ;
fprintf('All tests in %s.m passed.\n', this_script_name) ;
