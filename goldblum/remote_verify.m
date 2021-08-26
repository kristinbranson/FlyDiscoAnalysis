function remote_verify(source_user, source_host, source_path, dest_path)
    % record the start time
    tic_id = tic() ;

    % print an informative message
    fprintf('Verifying that contents of\n%s@%s:%s\nare present in\n%s\n... ', source_user, source_host, source_path, dest_path) ;
        
    % call helper
    % All those zeros are the numbers of different kinds of things that have been verified so far
    spinner = spinner_object() ;
    [n_files_verified, n_dirs_verified, n_file_bytes_verified] = ...
        remote_verify_helper(source_user, source_host, source_path, dest_path, 0, 0, 0, spinner) ;  %#ok<ASGLU>  % this will throw if there's a verification failure
    spinner.stop() ;    
    
    % print the number of files etc verified
    %fprintf("%d files verified\n", n_files_verified) ;
    %fprintf("%d folders verified\n", n_dirs_verified) ;
    %fprintf("%d file bytes verified\n", n_file_bytes_verified) ;
    fprintf("Success: All files and folders in source are present in destination.\n")

    % print the elapsed time
    elapsed_time = toc(tic_id) ;
    fprintf("Elapsed time: %0.1f seconds\n", elapsed_time) ;
end
