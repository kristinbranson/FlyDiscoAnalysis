function relative_path_from_synched_experiment_index = ...
        remote_sync_verify_and_delete_experiment_folders(source_user_name, ...
                                                         source_host_name, ...
                                                         source_root_absolute_path, ...
                                                         dest_root_absolute_path, ...
                                                         to_process_folder_name)
    % Make sure the remote folder exists, return if not
    does_folder_exist = does_remote_file_exist(source_user_name, source_host_name, source_root_absolute_path) ;
    if ~does_folder_exist ,
        fprintf('Folder %s does not exist on host %s, so not searching for experiment folders in it.\n', source_root_absolute_path, source_host_name) ;
        relative_path_from_synched_experiment_index = cell(0,1) ;
        return
    end
    
    % record the start time
    tic_id = tic() ;
        
    % print an informative message
    %fprintf('Searching for experiment folders in\n  %s@%s:%s...', source_user_name, source_host_name, source_root_absolute_path) ;
    [relative_path_from_experiment_folder_index, is_aborted_from_experiment_folder_index] = ...
        find_remote_experiment_folders(source_user_name, source_host_name, source_root_absolute_path, to_process_folder_name) ;
    experiment_folder_count = length(relative_path_from_experiment_folder_index) ;
    %fprintf("%d experiment folders found.\n" , experiment_folder_count) ;
    
    % sort the aborted from unaborted experiments
    relative_path_from_aborted_experiment_folder_index = relative_path_from_experiment_folder_index(is_aborted_from_experiment_folder_index) ;
    relative_path_from_unaborted_experiment_folder_index = relative_path_from_experiment_folder_index(~is_aborted_from_experiment_folder_index) ;    

    % print an informative message
    aborted_experiment_folder_count = length(relative_path_from_aborted_experiment_folder_index) ;
    if aborted_experiment_folder_count==1 ,
        fprintf('Deleting %d ABORTED experiment folder from\n  %s@%s:%s\n  ...\n', ...
                aborted_experiment_folder_count, source_user_name, source_host_name, source_root_absolute_path) ;
    else
        fprintf('Deleting %d ABORTED experiment folders from\n  %s@%s:%s\n  ...\n', ...
                aborted_experiment_folder_count, source_user_name, source_host_name, source_root_absolute_path) ;
    end
    
    % Delete each experiment folder in turn
    did_delete_from_aborted_experiment_folder_index = false(experiment_folder_count, 1) ;
    for i = 1 : aborted_experiment_folder_count ,
        experiment_folder_relative_path = relative_path_from_aborted_experiment_folder_index{i} ;
        source_folder_absolute_path = fullfile(source_root_absolute_path, experiment_folder_relative_path) ;
        try
            delete_remote_folder(source_user_name, source_host_name, source_folder_absolute_path) ;
            did_delete_from_aborted_experiment_folder_index(i) = true ;
        catch me ,
            fprintf('There was a problem during the deleting of ABORTED source experiment folder\n  %s\nThe problem was:\n%s\n', ... 
                    source_folder_absolute_path, ...
                    me.getReport()) ;    
        end            
    end

    % print the number of ABORTED experiment folders deleted
    deleted_aborted_experiment_folder_count = sum(double(did_delete_from_aborted_experiment_folder_index)) ;
    delete_error_count = aborted_experiment_folder_count - deleted_aborted_experiment_folder_count ;
    fprintf("Of %d ABORTED experiment folders:\n", aborted_experiment_folder_count) ;
    fprintf("  %d deleted\n", deleted_aborted_experiment_folder_count) ;
    fprintf("  %d failed to delete\n", delete_error_count) ;
    
    % print an informative message
    unaborted_experiment_folder_count = length(relative_path_from_unaborted_experiment_folder_index) ;
    fprintf('Synching %d experiment folders from\n  %s@%s:%s\n  into\n  %s\n  ...\n', ...
            unaborted_experiment_folder_count, source_user_name, source_host_name, source_root_absolute_path, dest_root_absolute_path) ;

    % Sync each experiment folder in turn
    did_synch_from_unaborted_experiment_folder_index = false(unaborted_experiment_folder_count, 1) ;
    for i = 1 : unaborted_experiment_folder_count ,
        experiment_folder_relative_path = relative_path_from_unaborted_experiment_folder_index{i} ;
        source_folder_absolute_path = fullfile(source_root_absolute_path, experiment_folder_relative_path) ;
        dest_folder_absolute_path = fullfile(dest_root_absolute_path, experiment_folder_relative_path) ;
        try
            remote_sync_and_verify(source_user_name, ...
                                   source_host_name, ...
                                   source_folder_absolute_path, ...
                                   dest_folder_absolute_path) ;
            did_synch_from_unaborted_experiment_folder_index(i) = true ;
        catch me ,
            fprintf('There was a problem during the synch of source experiment folder\n  %s\nThe problem was:\n%s\n', ... 
                    source_folder_absolute_path, ...
                    me.getReport()) ;    
        end            
    end

    % print the number of experiment folders copied
    synched_experiment_folder_count = sum(double(did_synch_from_unaborted_experiment_folder_index)) ;
    synch_error_count = unaborted_experiment_folder_count - synched_experiment_folder_count ;
    fprintf("Of %d unaborted experiment folders:\n", unaborted_experiment_folder_count) ;
    fprintf("  %d synched and verified\n", synched_experiment_folder_count) ;
    fprintf("  %d failed to synch or verify\n", synch_error_count) ;
    
    % Delete each synched experiment folder in turn
    relative_path_from_synched_experiment_index = relative_path_from_unaborted_experiment_folder_index(did_synch_from_unaborted_experiment_folder_index) ;
    synched_experiment_count = length(relative_path_from_synched_experiment_index) ;
    did_delete_from_synched_experiment_index = false(synched_experiment_count, 1) ;
    for i = 1 : unaborted_experiment_folder_count ,
        if did_synch_from_unaborted_experiment_folder_index(i) ,
            experiment_folder_relative_path = relative_path_from_unaborted_experiment_folder_index{i} ;
            source_folder_absolute_path = fullfile(source_root_absolute_path, experiment_folder_relative_path) ;
            try
                delete_remote_folder(source_user_name, source_host_name, source_folder_absolute_path) ;
                did_delete_from_synched_experiment_index(i) = true ;
            catch me ,
                fprintf('There was a problem during the post-synch deletion of source experiment folder\n  %s\nThe problem was:\n%s\n', ...
                    source_folder_absolute_path, ...
                    me.getReport()) ;
            end
        end
    end
    
    % print the number of experiment folders copied
    deleted_experiment_folder_count = sum(double(did_delete_from_synched_experiment_index)) ;
    delete_error_count = synched_experiment_folder_count - deleted_experiment_folder_count ;
    fprintf("Of %d synched experiment folders:\n", synched_experiment_folder_count) ;
    fprintf("  %d deleted\n", deleted_experiment_folder_count) ;
    fprintf("  %d failed to delete\n", delete_error_count) ;

    % print the elapsed time
    elapsed_time = toc(tic_id) ;
    fprintf("Total elapsed time: %0.1f seconds\n", elapsed_time) ;
    
%     % throw an error if there were any failures
%     if synch_error_count > 0 || delete_error_count > 0 ,
%         error("There was at least one failure during the synching of unaborted experiment folders from the remote host") ;
%     end
    
    % Return the synched experiments, whether or not they were deleted
    relative_path_from_synched_experiment_index = relative_path_from_unaborted_experiment_folder_index(did_synch_from_unaborted_experiment_folder_index) ;
end
