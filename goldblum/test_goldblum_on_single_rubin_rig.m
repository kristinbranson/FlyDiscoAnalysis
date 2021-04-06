do_use_bqueue = true ;
do_actually_submit_jobs = true ;
do_run_analysis_in_debug_mode = true ;

% Specify the "per-lab" configuration here
per_lab_configuration = struct() ;
per_lab_configuration.lab_head_last_name = 'scicompsoft' ;
per_lab_configuration.host_name_from_rig_index = {'flybowl-ww1.hhmi.org'} ;
per_lab_configuration.rig_user_name_from_rig_index = {'labadmin'} ;
per_lab_configuration.data_folder_path_from_rig_index = {'/cygdrive/h/flydisco_data'} ;
per_lab_configuration.destination_folder = '/groups/branson/bransonlab/taylora/flydisco/goldblum/goldblum-test-destination-folder' ;    
per_lab_configuration.settings_folder_path = '/groups/branson/bransonlab/taylora/flydisco/goldblum/FlyDiscoAnalysis/settings' ;
per_lab_configuration.does_use_per_user_folders = true ;

this_script_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_script_path) ;
analysis_test_template_folder_path = fullfile(this_folder_path, 'analysis-test-template') ;

% Delete the destination folder
destination_folder_path = fullfile(this_folder_path, 'goldblum-test-destination-folder') ;
if exist(destination_folder_path, 'file') ,
    return_code = system_with_error_handling(sprintf('rm -rf %s', destination_folder_path)) ;
end

% Recopy the analysis test folder from the template
reset_analysis_test_folder() ;

% Copy it to the rig computer
remote_host_name = 'flybowl-ww1.hhmi.org' ;
remote_user_folder_path = '/cygdrive/h/flydisco_data/scicompsoft/taylora' ;
remote_user_name = 'labadmin' ;
command_line = sprintf('scp -r %s %s@%s:%s', analysis_test_template_folder_path, remote_user_name, remote_host_name, remote_user_folder_path) ;
system_with_error_handling(command_line) ;

% Run goldblum
goldblum(do_use_bqueue, do_actually_submit_jobs, do_run_analysis_in_debug_mode, per_lab_configuration) ;        

% Check that the expected files are present on dm11
analysis_test_folder_destination_path = fullfile(destination_folder_path, 'taylora', 'analysis-test-folder') ;
local_verify(analysis_test_template_folder_path, analysis_test_folder_destination_path) ;

% Check that some of the expected outputs were generated
test_file_names = {'perframe' 'scores_AttemptedCopulation.mat' 'scoresBackup.mat' 'registered_trx.mat' 'wingtracking_results.mat'} ;
for i = 1 : length(test_file_names) ,
    test_file_name = test_file_names{i} ;
    test_file_path = ...
        fullfile(analysis_test_folder_destination_path, ...
                 '2020-01-07', ...
                 'SS36564_20XUAS_CsChrimson_mVenus_attP18_flyBowlMing_20200227_Continuous_2min_5int_20200107_20200229T132141', ...
                 test_file_name) ;
    if ~exist(test_file_path, 'file') ,
        error('No output file at %s', test_file_path) ;
    end
end    

% Check that the rig user folder is empty now
entry_names = ...
    list_remote_dir(remote_user_name, remote_host_name, remote_user_folder_path) ;
if ~isempty(entry_names) ,
    error('Remote user folder %s:%s is not empty', remote_host_name, remote_user_folder_path) ;
end

% Run goldblum again, make sure nothing has changed
goldblum(do_use_bqueue, do_actually_submit_jobs, do_run_analysis_in_debug_mode, per_lab_configuration) ;        

% Check that the expected files are present on dm11
analysis_test_folder_destination_path = fullfile(destination_folder_path, 'taylora', 'analysis-test-folder') ;
local_verify(analysis_test_template_folder_path, analysis_test_folder_destination_path) ;

% Check that some of the expected outputs were generated
test_file_names = {'perframe' 'scores_AttemptedCopulation.mat' 'scoresBackup.mat' 'registered_trx.mat' 'wingtracking_results.mat'} ;
for i = 1 : length(test_file_names) ,
    test_file_name = test_file_names{i} ;
    test_file_path = ...
        fullfile(analysis_test_folder_destination_path, ...
                 '2020-01-07', ...
                 'SS36564_20XUAS_CsChrimson_mVenus_attP18_flyBowlMing_20200227_Continuous_2min_5int_20200107_20200229T132141', ...
                 test_file_name) ;
    if ~exist(test_file_path, 'file') ,
        error('No output file at %s', test_file_path) ;
    end
end    

% Check that the rig user folder is empty now
entry_names = ...
    list_remote_dir(remote_user_name, remote_host_name, remote_user_folder_path) ;
if ~isempty(entry_names) ,
    error('Remote user folder %s:%s is not empty', remote_host_name, remote_user_folder_path) ;
end

% If get here, all is well
[~, this_script_name] = fileparts(this_script_path) ;
fprintf('All tests in %s.m passed.\n', this_script_name) ;
