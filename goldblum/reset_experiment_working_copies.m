function reset_experiment_working_copies(working_copy_experiments_folder_path, read_only_experiments_folder_path)
    if exist(working_copy_experiments_folder_path, 'file') ,
        system_from_list_with_error_handling({'rm', '-rf', working_copy_experiments_folder_path}) ;
    end
    system_from_list_with_error_handling( ...
        {'cp', '--no-preserve=mode', '-R', '-T', read_only_experiments_folder_path, working_copy_experiments_folder_path} ) ;
end
