% Where does this script live?
this_script_path = mfilename('fullpath') ;
fly_disco_analysis_test_folder_path = fileparts(this_script_path) ;
fly_disco_analysis_folder_path = fileparts(fly_disco_analysis_test_folder_path) ;
fly_disco_folder_path = fileparts(fly_disco_analysis_folder_path) ;

% Set some params
settings_folder_path = fullfile(fly_disco_analysis_folder_path, 'settings-internal') ;  % for now, want to use internal settings
analysis_protocol = '20240402_flybubble_LED_LPC1_CsChr' ;
do_reset_working_experiments_folder = true ;
submit_host_name = if_not_a_submit_host('submit.int.janelia.org') ;

% Set the name of the read-only folder and the working folder
read_only_experiments_folder_path = fullfile(fly_disco_folder_path, 'example-experiments', 'stimoff-issue-read-only') ;
working_experiments_folder_path = fullfile(fly_disco_folder_path, 'example-experiments', 'stimoff-issue') ;

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
experiment_count = numel(folder_path_from_experiment_index) ;
for experiment_index = 1 : experiment_count ,
  expdir = folder_path_from_experiment_index{experiment_index} ;
  FlyDiscoComputePerFrameStats(expdir,'settingsdir',settings_folder_path,'analysis_protocol',analysis_protocol,'dorecompute',false,'docomputehists',true,...
                               'debugplot',false);
  FlyDiscoPlotPerFrameStats(expdir,'settingsdir',settings_folder_path,'analysis_protocol',analysis_protocol,'debug',false,'makestimvideos',1, ...
                            'plotstim',true,'plothist',true,'plotflies',true,'plotstimtrajs',true,'plottimeseries',1);
end
