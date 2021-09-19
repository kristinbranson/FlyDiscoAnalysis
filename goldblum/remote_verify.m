function remote_verify(source_user, source_host, source_path, dest_path, be_verbose)
    % Deal with arguments
    if ~exist('be_verbose', 'var') || isempty(be_verbose) ,
        be_verbose = false ;
    end
    
    % record the start time
    tic_id = tic() ;

    % print an informative message
    if be_verbose ,
        fprintf('Verifying that contents of\n%s@%s:%s\nare present in\n%s\n... ', source_user, source_host, source_path, dest_path) ;
    end
    
    % call helper
    % All those zeros are the numbers of different kinds of things that have been verified so far
    if be_verbose ,
        spinner = spinner_object() ;
    else
        spinner = spinner_object('mute') ;
    end        
    [n_files_verified, n_dirs_verified, n_file_bytes_verified] = ...
        remote_verify_helper(source_user, source_host, source_path, dest_path, 0, 0, 0, spinner) ;  %#ok<ASGLU>  
            % this will throw if there's a verification failure
    spinner.stop() ;    
    
    % print the number of files etc verified
    %fprintf("%d files verified\n", n_files_verified) ;
    %fprintf("%d folders verified\n", n_dirs_verified) ;
    %fprintf("%d file bytes verified\n", n_file_bytes_verified) ;
    if be_verbose ,
        fprintf("Success: All files and folders in source are present in destination.\n")
    end

    % print the elapsed time
    elapsed_time = toc(tic_id) ;
    if be_verbose ,
        fprintf("Elapsed time: %0.1f seconds\n", elapsed_time) ;
    end
end
