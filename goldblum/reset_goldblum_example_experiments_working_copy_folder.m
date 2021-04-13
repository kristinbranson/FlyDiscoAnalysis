function working_copy_example_experiment_folder_path = reset_goldblum_example_experiments_working_copy_folder()
    this_script_path = mfilename('fullpath') ;
    this_folder_path = fileparts(this_script_path) ;
    read_only_example_experiment_folder_path = fullfile(this_folder_path, 'example-experiments-read-only') ;
    working_copy_example_experiment_folder_path = fullfile(this_folder_path, 'example-experiments-working-copy') ;
    if exist(working_copy_example_experiment_folder_path, 'file') ,
        rmdir(working_copy_example_experiment_folder_path, 's') ;
    end
    did_succeed = mkdir(working_copy_example_experiment_folder_path) ;
    if ~did_succeed ,
        error('Unable to create folder %s', working_copy_example_experiment_folder_path) ;
    end
    did_succeed = copyfile(read_only_example_experiment_folder_path, working_copy_example_experiment_folder_path) ;
    if ~did_succeed ,
        error('Unable to copy %s to %s', read_only_example_experiment_folder_path, working_copy_example_experiment_folder_path) ;
    end
end
