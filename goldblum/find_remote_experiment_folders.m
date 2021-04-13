function experiment_folder_relative_path_list = find_remote_experiment_folders(user_name, host_name, path, to_process_folder_name)
    % record the start time
    start_time = tic() ;

    % print an informative message
    fprintf("Looking for experiment folders within %s on host %s as user %s... ", path, host_name, user_name)
    
    % call helper
    % All those zeros are the numbers of different kinds of things that have been verified so far
    spinner = spinner_object() ;
    experiment_folder_relative_path_list = find_remote_experiment_folders_helper(user_name, host_name, '', path, to_process_folder_name, spinner) ;
    spinner.stop() ;
   
    % print the number of experiment folders found
    fprintf("%d experiment folders found\n" , length(experiment_folder_relative_path_list))

    % print the elapsed time
    elapsed_time = toc(start_time) ;
    fprintf("Elapsed time: %0.1f seconds\n", elapsed_time) ;
end
