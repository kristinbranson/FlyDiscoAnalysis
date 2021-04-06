function remote_sync(source_user, source_host, source_path, dest_path)
    % record the start time
    tic_id = tic() ;

    % if dest dir doesn't exist, create it
    if exist(dest_path, 'file') ,
        % make sure it's a dir
        if ~exist(dest_path, 'dir') ,
            error('Destination %s exists, but is a file, not a directory', dest_path)
        end
    else
        ensure_folder_exists(dest_path) ;
    end
        
    % print an informative message
    fprintf('Copying contents of\n%s@%s:%s\ninto\n%s\n... ', source_user, source_host, source_path, dest_path) ;
    
    % call helper
    % All those zeros are the numbers of different kinds of failures so far
    spinner = spinner_object() ;
    [n_copied, n_failed, n_dir_failed, n_verified, n_dir_failed_to_list, time_spent_copying] = ...
        remote_sync_helper(source_user, source_host, source_path, dest_path, 0, 0, 0, 0, 0, 0.0, spinner) ;
    spinner.stop() ;
    
    % print the number of files copied
    fprintf("%d files copied\n", n_copied) ;
    fprintf("%d files verified\n", n_verified) ;
    fprintf("%d files copied or verified\n", (n_copied+n_verified)) ;
    fprintf("%d files failed to copy\n", n_failed) ;
    fprintf("%d directories failed to copy\n", n_dir_failed) ;
    fprintf("%d directories failed to list\n", n_dir_failed_to_list) ;

    % print the elapsed time
    elapsed_time = toc(tic_id) ;
    fprintf("Elapsed time: %0.1f seconds\n", elapsed_time) ;
    fprintf("Time spent copying: %0.1f seconds\n", time_spent_copying) ;

    % throw an error if there were any failures
    if n_failed + n_dir_failed + n_dir_failed_to_list > 0 ,
        error("There was at least one failure during the remote sync: %d file copies failed, %d directory creates failed, %d directories failed to list", ... 
              n_failed, n_dir_failed, n_dir_failed_to_list) ;
    end
end
