function experiment_folder_path_list = find_experiment_folders_relative(source_path)
    % record the start time
    start_time = tic() ;

    % Make sure the source_path exits
    if exist(source_path, 'file') ,
        if exist(source_path, 'dir') ,
          % All is well
        else
          error('%s is a file, not a folder', source_path) ;
        end
    else
        error('%s does not exist', source_path) ;
    end

    % print an informative message
    fprintf("Looking for experiment folders within %s... ", source_path)
    
    % call helper
    % All those zeros are the numbers of different kinds of things that have been verified so far
    spinner = spinner_object() ;
    experiment_folder_path_list = find_experiment_folders_relative_helper(source_path, '', spinner) ;
    spinner.stop() ;
   
    % print the number of files etc verified
    fprintf("%d experiment folders found\n" , length(experiment_folder_path_list))

    % print the elapsed time
    elapsed_time = toc(start_time) ;
    fprintf("Elapsed time: %0.1f seconds\n", elapsed_time) ;
end
