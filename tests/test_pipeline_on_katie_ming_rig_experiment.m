experiment_name = 'emptysplit_20xUAS-ChrimsonRmVenusattp18_flyBowlMing_nopause_lengthofpersis_2min_10int_20191218T093239_2' ;
reset_test_experiment_folder(experiment_name);

this_script_file_path = mfilename('fullpath') ;
this_script_folder_path = fileparts(this_script_file_path) ;
experiment_folder_path = fullfile(this_script_folder_path, 'analysis-test-folder', experiment_name) ;
settings_folder_path = fullfile(this_script_folder_path, 'settings-internal') ;

% Read the experiment metadata to determine the analysis_protoocol
metadata_file_path = determine_metadata_file_path(experiment_folder_path) ;
analysis_protocol_folder_name = analysis_protocol_from_metadata_file(metadata_file_path, settings_folder_path) ;
%fprintf('Analysis protocol is: %s\n', analysis_protocol_folder_name) ;

% Build up the parameters cell array
params = ...
    {'settingsdir',settings_folder_path,...
     'analysis_protocol',analysis_protocol_folder_name, ...
     'forcecompute',false,...
     'doautomaticchecksincoming',true,...
     'doflytracking',true,...
     'doregistration',true,...
     'doledonoffdetection',true,...
     'dotrackwings',true,...
     'dosexclassification',true,...
     'docomputeperframefeatures',true,...
     'docomputehoghofperframefeatures',false,...
     'dojaabadetect',true,...
     'docomputeperframestats',false,...
     'doplotperframestats',false,...
     'domakectraxresultsmovie',false,...
     'doextradiagnostics',false,...
     'doanalysisprotocol',isunix(),...
     'doautomaticcheckscomplete',false } ;

% Call the function to do the real work
FlyDiscoPipeline(experiment_folder_path, params{:}) ;
