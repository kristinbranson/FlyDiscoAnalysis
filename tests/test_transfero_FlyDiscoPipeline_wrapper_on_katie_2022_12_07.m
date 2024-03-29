% This isn't a proper test b/c it doesn't check whether anything worked

user_name_for_configuration_purposes = 'rubinlab' ;                                
analysis_parameters = {'do_try', true} ;
% analysis_parameters = ...
%          {'doautomaticchecksincoming',true,...
%           'doflytracking',true, ...
%           'doregistration',true,...
%           'doledonoffdetection',true,...
%           'dosexclassification',true,...
%           'dotrackwings',false,...
%           'docomputeperframefeatures',true,...
%           'docomputehoghofperframefeatures',false,...
%           'dojaabadetect',true,...
%           'docomputeperframestats',false,...
%           'doplotperframestats',false,...
%           'domakectraxresultsmovie',true,...
%           'doextradiagnostics',false,...
%           'doanalysisprotocol',true,...
%           'doautomaticcheckscomplete',false};

% Where does this script live?
this_script_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_script_path) ;
fly_disco_analysis_folder_path = fileparts(this_folder_path) ;
fly_disco_folder_path = fileparts(fly_disco_analysis_folder_path) ;
%settings_folder_path = fullfile(fly_disco_analysis_folder_path, 'settings-internal') ;
read_only_experiments_folder_path = fullfile(fly_disco_folder_path, 'example-experiments', 'single-katie-experiment-2022-12-07-read-only') ;
working_experiments_folder_path = fullfile(fly_disco_folder_path, 'example-experiments', 'single-katie-experiment-2022-12-07') ;

% Delete the destination folder
if exist(working_experiments_folder_path, 'file') ,
    return_code = system_from_list_with_error_handling({'rm', '-rf', working_experiments_folder_path}) ;
end

% Recopy the test folder from the template
fprintf('Resetting working experiments folder...\n') ;
reset_experiment_working_copies(working_experiments_folder_path, read_only_experiments_folder_path) ;

% Find the experiments
folder_path_from_experiment_index = find_experiment_folders(working_experiments_folder_path) ;

% Run the script under test
experiment_count = length(folder_path_from_experiment_index) ;
for experiment_index = 1 : experiment_count ,
    experiment_folder_path = folder_path_from_experiment_index{experiment_index} ;
    fprintf('Running transfero_FlyDiscoPipeline_wrapper on experiment index %d...\n', experiment_index) ;
    transfero_FlyDiscoPipeline_wrapper(experiment_folder_path, user_name_for_configuration_purposes, analysis_parameters{:}) ;
end
