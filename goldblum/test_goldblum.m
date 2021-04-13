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
destination_folder_path = fullfile(this_folder_path, 'goldblum-test-destination-folder') ;
settings_folder_path = fullfile(fly_disco_analysis_folder_path, 'settings') ;
per_lab_configuration = struct() ;
per_lab_configuration.lab_head_last_name = lab_head_last_name ;
per_lab_configuration.host_name_from_rig_index = {rig_host_name} ;
per_lab_configuration.rig_user_name_from_rig_index = {rig_user_name} ;
per_lab_configuration.data_folder_path_from_rig_index = {rig_data_folder_path} ;
per_lab_configuration.destination_folder = destination_folder_path ;    
per_lab_configuration.settings_folder_path = settings_folder_path ;
%per_lab_configuration.does_use_per_user_folders = true ;

% Delete the destination folder
if exist(destination_folder_path, 'file') ,
    return_code = system_with_error_handling(sprintf('rm -rf %s', destination_folder_path)) ;
end

% Recopy the analysis test folder from the template
fprintf('Resetting analysis test folder...\n') ;
example_experiments_folder_path = reset_goldblum_example_experiments_working_copy_folder() ;

% Copy it to the rig computer (or just direct to the destination folder)
fprintf('Transfering data to the rig computer (or to the destination path if not transfering data from rig)...\n') ;
rig_lab_data_folder_path = fullfile(rig_data_folder_path, lab_head_last_name) ;
example_experiments_folder_destination_path = fullfile(destination_folder_path) ;
if do_transfer_data_from_rigs ,
    command_line = sprintf('scp -B -r %s/* %s@%s:%s', example_experiments_folder_path, rig_user_name, rig_host_name, rig_lab_data_folder_path) ; %#ok<UNRCH>
else
    ensure_folder_exists(example_experiments_folder_destination_path) ;  %#ok<UNRCH>
    command_line = sprintf('cp -R %s/* %s', example_experiments_folder_path, example_experiments_folder_destination_path) ;
end
system_with_error_handling(command_line) ;

% Run goldblum
fprintf('Running goldblum...\n') ;
goldblum(do_transfer_data_from_rigs, do_run_analysis, do_use_bqueue, do_actually_submit_jobs, [], per_lab_configuration) ;        

% Check that the expected files are present on dm11
local_verify(example_experiments_folder_path, example_experiments_folder_destination_path) ;

% Check that some of the expected outputs were generated
test_file_names = {'perframe' 'scores_AttemptedCopulation.mat' 'scoresBackup.mat' 'registered_trx.mat' 'wingtracking_results.mat'} ;
for i = 1 : length(test_file_names) ,
    test_file_name = test_file_names{i} ;
    test_file_path = ...
        fullfile(example_experiments_folder_destination_path, ...
                 'SS36564_20XUAS_CsChrimson_mVenus_attP18_flyBowlMing_20200227_Continuous_2min_5int_20200107_20200229T132141', ...
                 test_file_name) ;
    if ~exist(test_file_path, 'file') ,
        error('No output file at %s', test_file_path) ;
    end
end    

% Check that the rig user folder is empty now
entry_names = ...
    list_remote_dir(rig_user_name, rig_host_name, rig_lab_data_folder_path) ;
if ~isempty(entry_names) ,
    error('Rig lab data folder %s:%s is not empty', rig_host_name, rig_lab_data_folder_path) ;
end

% Run goldblum again, make sure nothing has changed
goldblum(do_transfer_data_from_rigs, do_use_bqueue, do_actually_submit_jobs, [], per_lab_configuration) ;        

% Check that the expected files are present on dm11
example_experiments_folder_destination_path = fullfile(destination_folder_path, 'taylora', 'analysis-test-folder') ;
local_verify(example_experiments_folder_path, example_experiments_folder_destination_path) ;

% Check that some of the expected outputs were generated
test_file_names = {'perframe' 'scores_AttemptedCopulation.mat' 'scoresBackup.mat' 'registered_trx.mat' 'wingtracking_results.mat'} ;
for i = 1 : length(test_file_names) ,
    test_file_name = test_file_names{i} ;
    test_file_path = ...
        fullfile(example_experiments_folder_destination_path, ...
                 'SS36564_20XUAS_CsChrimson_mVenus_attP18_flyBowlMing_20200227_Continuous_2min_5int_20200107_20200229T132141', ...
                 test_file_name) ;
    if ~exist(test_file_path, 'file') ,
        error('No output file at %s', test_file_path) ;
    end
end    

% Check that the rig user folder is empty now
entry_names = ...
    list_remote_dir(rig_user_name, rig_host_name, rig_lab_data_folder_path) ;
if ~isempty(entry_names) ,
    error('Remote user folder %s:%s is not empty', rig_host_name, rig_lab_data_folder_path) ;
end

% If get here, all is well
[~, this_script_name] = fileparts(this_script_path) ;
fprintf('All tests in %s.m passed.\n', this_script_name) ;
