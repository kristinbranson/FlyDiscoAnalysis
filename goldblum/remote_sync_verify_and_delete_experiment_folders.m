function result = remote_sync_verify_and_delete_experiment_folders(source_user_name, ...
                                                                   source_host_name, ...
                                                                   source_root_absolute_path, ...
                                                                   dest_root_absolute_path, ...
                                                                   to_process_folder_name)
    % Make sure the remote folder exists, return if not
    does_folder_exist = does_remote_file_exist(source_user_name, source_host_name, source_root_absolute_path) ;
    if ~does_folder_exist ,
        fprintf('Folder %s does not exist on host %s, so not searching for experiment folders in it.\n', source_root_absolute_path, source_host_name) ;
        result = cell(0,1) ;
        return
    end
    
    % record the start time
    tic_id = tic() ;
        
    % print an informative message
    fprintf('Searching for experiment folders in\n%s@%s:%s\n... ', source_user_name, source_host_name, source_root_absolute_path) ;
    relative_path_from_experiment_folder_index = ...
        find_remote_experiment_folders(source_user_name, source_host_name, source_root_absolute_path, to_process_folder_name) ;
    experiment_folder_count = length(relative_path_from_experiment_folder_index) ;
    
    % print an informative message
    fprintf('Copying %d experiment folders from\n%s@%s:%s\ninto\n%s\n... ', ...
            experiment_folder_count, source_user_name, source_host_name, source_root_absolute_path, dest_root_absolute_path) ;

    % Sync each experiment folder in turn
    did_synch_from_experiment_folder_index = false(1, experiment_folder_count) ;
    for i = 1 : experiment_folder_count ,
        experiment_folder_relative_path = relative_path_from_experiment_folder_index{i} ;
        source_folder_absolute_path = fullfile(source_root_absolute_path, experiment_folder_relative_path) ;
        dest_folder_absolute_path = fullfile(dest_root_absolute_path, experiment_folder_relative_path) ;
        try
            remote_sync_and_verify_and_delete(source_user_name, ...
                                              source_host_name, ...
                                              source_folder_absolute_path, ...
                                              dest_folder_absolute_path) ;
            did_synch_from_experiment_folder_index(i) = true ;
        catch me ,
            fprintf('There was a problem during the synch of source experiment folder\n  %s\nThe problem was:\n%s\n', ... 
                    source_folder_absolute_path, ...
                    me.getReport()) ;    
        end            
    end
    
    % print the number of experiment folders copied
    synched_experiment_folder_count = sum(double(did_synch_from_experiment_folder_index)) ;
    error_count = experiment_folder_count - synched_experiment_folder_count ;
    fprintf("%d experiment folders synched, verified, and deleted\n", synched_experiment_folder_count) ;
    fprintf("%d experiment folders failed to synch, verify, or delete\n", error_count) ;

    % print the elapsed time
    elapsed_time = toc(tic_id) ;
    fprintf("Elapsed time: %0.1f seconds\n", elapsed_time) ;
    
    % throw an error if there were any failures
    if error_count > 0 ,
        error("There was at least one failure during the synching of experiment folders from the remote host") ;
    end
    
    % Return the synched experiments
    result = relative_path_from_experiment_folder_index(did_synch_from_experiment_folder_index) ;
end
