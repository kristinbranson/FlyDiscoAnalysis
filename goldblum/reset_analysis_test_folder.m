function reset_analysis_test_folder()
    this_script_path = mfilename('fullpath') ;
    this_folder_path = fileparts(this_script_path) ;
    template_test_experiment_folder_path = fullfile(this_folder_path, 'analysis-test-template') ;
    test_experiment_folder_path = fullfile(this_folder_path, 'analysis-test-folder') ;
    if exist(test_experiment_folder_path, 'file') ,
        rmdir(test_experiment_folder_path, 's') ;
    end
    did_succeed = mkdir(test_experiment_folder_path) ;
    if ~did_succeed ,
        error('Unable to create folder %s', test_experiment_folder_path) ;
    end
    did_succeed = copyfile(template_test_experiment_folder_path, test_experiment_folder_path) ;
    if ~did_succeed ,
        error('Unable to copy %s to %s', template_test_experiment_folder_path, test_experiment_folder_path) ;
    end
end
