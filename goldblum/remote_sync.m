function remote_sync(source_user, source_host, source_path, dest_path, be_verbose)
    % Deal with arguments
    if ~exist('be_verbose', 'var') || isempty(be_verbose) ,
        be_verbose = false ;
    end
    
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
    if be_verbose ,
        fprintf('Copying contents of\n%s@%s:%s\ninto\n%s\n... ', source_user, source_host, source_path, dest_path) ;
    end
    
    % call helper
    % All those zeros are the numbers of different kinds of failures so far
    if be_verbose ,
        spinner = spinner_object() ;
    else
        spinner = spinner_object('mute') ;
    end
    [n_copied, n_failed, n_dir_failed, n_verified, n_dir_failed_to_list, time_spent_copying] = ...
        remote_sync_helper(source_user, source_host, source_path, dest_path, 0, 0, 0, 0, 0, 0.0, spinner) ;  %#ok<ASGLU>
    spinner.stop() ;
    
    % print the number of files copied
    if n_failed==0 && n_dir_failed==0 && n_dir_failed_to_list==0 ,
        if be_verbose ,
            fprintf('Successfully copied contents of\n%s@%s:%s\ninto\n%s\n', source_user, source_host, source_path, dest_path) ;
        end

        % print the elapsed time
        elapsed_time = toc(tic_id) ;
        if be_verbose ,
            fprintf("Elapsed time: %0.1f seconds\n", elapsed_time) ;
            fprintf("Time spent copying: %0.1f seconds\n", time_spent_copying) ;
        end
    else
        % throw an error if there were any failures
        error('remote_sync:did_fail', ...
              'There was at least one failure during the remote sync: %d file copies failed, %d directory creates failed, %d directories failed to list', ... 
              n_failed, n_dir_failed, n_dir_failed_to_list) ;
    end
end
