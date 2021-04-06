reset_test_experiment_folder();

this_script_file_path = mfilename('fullpath') ;
this_script_folder_path = fileparts(this_script_file_path) ;
experiment_folder_path = ...
    fullfile(this_script_folder_path, ...
             'analysis-test-folder', ...
             'locomotionGtACR1_24_EXT_VGLUT-GAL4_RigA_20210226T095136') ;
analysis_protocol_folder_path = fullfile(this_script_folder_path, 'FlyDiscoAnalysis/settings/20210228_flybubble_dickson_locomotionGtACR1') ;

do_run_analysis_in_debug_mode = false ;
fly_disco_analysis_pipeline_wrapper(experiment_folder_path, analysis_protocol_folder_path, do_run_analysis_in_debug_mode) ;
