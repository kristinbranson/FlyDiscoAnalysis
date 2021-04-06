function local_verify(source_path, dest_path)
    % record the start time
    tic_id = tic() ;

    % print an informative message
    fprintf("Verifying that contents of\n%s\nare also present in\n%s\n... ", source_path, dest_path) ;
    
    % call helper
    % All those zeros are the numbers of different kinds of things that have been verified so far
    spinner = spinner_object() ;
    [n_files_verified, n_dirs_verified, n_file_bytes_verified] = ...
        local_verify_helper(source_path, dest_path, 0, 0, 0, spinner) ;
    spinner.stop() ;    
    
    % print the number of files etc verified
    fprintf("%d files verified\n", n_files_verified) ;
    fprintf("%d folders verified\n", n_dirs_verified) ;
    fprintf("%d file bytes verified\n", n_file_bytes_verified) ;
    fprintf("Success: All files and folders in source are present in destination.\n")

    % print the elapsed time
    elapsed_time = toc(tic_id) ;
    fprintf("Elapsed time: %0.1f seconds\n", elapsed_time) ;
end
