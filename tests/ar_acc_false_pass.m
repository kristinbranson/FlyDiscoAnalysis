% Where does this script live?
this_script_path = mfilename('fullpath') ;
fly_disco_analysis_test_folder_path = fileparts(this_script_path) ;
fly_disco_analysis_folder_path = fileparts(fly_disco_analysis_test_folder_path) ;
fly_disco_folder_path = fileparts(fly_disco_analysis_folder_path) ;

cluster_billing_account_name = [] ;
user_name_for_configuration_purposes = 'bransonlab' ;
%analysis_parameters = cell(1,0) ;
analysis_parameters = ...
         {'docomputeperframefeatures',false, ...
          'doapt',false, ...
          'domakectraxresultsmovie',false } ;  % don't need these stages right now
% analysis_parameters = ...
%          { 'domakectraxresultsmovie', true } ;
settings_folder_path = fullfile(fly_disco_analysis_folder_path, 'settings-internal') ;  % for now, want to use internal settings
%settings_folder_path = '/groups/branson/bransonlab/taylora/flydisco/OtopalikFlyDiscoSettings/settings' ;  % at same commit as when prod run happened
do_use_bqueue = false ;
do_actually_submit_jobs = false ;
do_try = true ;
do_reset_working_experiments_folder = true ;
submit_host_name = if_not_a_submit_host('submit.int.janelia.org') ;
initial_optional_argument_list = { ...
  'settingsdir', settings_folder_path, ...
  'do_try', do_try, ...
  'debug', true } ;
optional_argument_list = horzcat(initial_optional_argument_list, analysis_parameters) ;

read_only_experiments_folder_path = fullfile(fly_disco_folder_path, 'example-experiments', 'ar-acc-false-pass-read-only') ;
working_experiments_folder_path = fullfile(fly_disco_folder_path, 'example-experiments', 'ar-acc-false-pass') ;

% Recopy the working folder from the read-only one
if do_reset_working_experiments_folder || ~exist(working_experiments_folder_path, 'dir')
  fprintf('Resetting working experiments folder...\n') ;
  tic_id = tic() ;
  reset_experiment_working_copies(working_experiments_folder_path, read_only_experiments_folder_path) ;
  elapsed_time = toc(tic_id) ;
  fprintf('Elapsed time: %g s\n', elapsed_time) ;
end

% Find the experiments
folder_path_from_experiment_index = find_experiment_folders(working_experiments_folder_path) ;

% Call the testing function to do the real work
run_transfero_FlyDiscoPipeline_wrapper_on_experiment_list(...
  folder_path_from_experiment_index, ...
  cluster_billing_account_name, ...
  user_name_for_configuration_purposes, ...
  do_use_bqueue, ...
  do_actually_submit_jobs, ...
  do_try, ...
  submit_host_name, ...
  optional_argument_list{:})

% All the ACC check files should exist, but for some experiments ACC check
% failure is the desired result.
does_acc_mat_file_exist_from_experiment_folder_index = cellfun(@does_acc_mat_file_exist, folder_path_from_experiment_index) ;
do_all_acc_mat_files_exist = all(did_succeed_from_experiment_folder_index) ;
if do_all_acc_mat_files_exist ,
  fprintf('Test failed---not all ACC mat files are present.\n') ;  
  return
end

% Check for success in the ACC output files
should_acc_pass_from_experiment_folder_index = cellfun(@determine_if_acc_checks_should_pass_on_experiment_folder, folder_path_from_experiment_index) ;
did_acc_pass_from_experiment_folder_index = cellfun(@did_acc_pass_for_experiment_folder, folder_path_from_experiment_index) ;
did_all_do_the_right_thing = all(did_acc_pass_from_experiment_folder_index==should_acc_pass_from_experiment_folder_index) ;
if did_all_do_the_right_thing ,
  fprintf('Test passed.\n') ;
else
  fprintf('Test failed---some files that should have passed ACC did not, or vice-versa.\n') ;  
end
