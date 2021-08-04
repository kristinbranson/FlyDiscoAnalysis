function reset_goldblum_example_experiments_working_copy_folder(example_experiment_folder_path, read_only_example_experiment_folder_path)
    if exist(example_experiment_folder_path, 'file') ,
        system_from_list_with_error_handling({'rm', '-rf', example_experiment_folder_path}) ;
    end
    system_from_list_with_error_handling( ...
        {'cp', '--no-preserve=mode', '-R', '-T', read_only_example_experiment_folder_path, example_experiment_folder_path} ) ;
end
