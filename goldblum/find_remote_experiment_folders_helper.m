function result = find_remote_experiment_folders_helper(user_name, host_name, parent_relative_path, root_absolute_path, to_process_folder_name,spinner)
    % Find the experiment folders on a remote host.  Returns realtive paths,
    % relative to root_absolute_path.
  
    % Get a list of all files and folders
    parent_absolute_path = fullfile(root_absolute_path, parent_relative_path) ;
    try
        [entries, ~, ~, is_entry_a_folder, ~] = ...
            list_remote_dir(user_name, host_name, parent_absolute_path) ;
    catch me ,
        % if we can't list the dir, warn but continue
        if isequal(me.identifier, 'list_remote_dir:failed') ,
            spinner.print("Warning: can't list path %s on host %s as user %s", parent_absolute_path, host_name, user_name) ;
            spinner.print("%s", me.getReport()) ;
            %n_dirs_failed_to_list = n_dirs_failed_to_list + 1 ;
            return
        else
            rethrow(me) ;
        end
    end
    spinner.spin() ;

    % Separate source file, folder names    
    file_names = entries(~is_entry_a_folder) ;
    folder_names = entries(is_entry_a_folder) ;    

    % If the parent_path is an experiment folder, we're done
    if is_experiment_folder_given_contents(file_names) ,
        if isequal(parent_relative_path, to_process_folder_name) ,
            spinner.print("Warning: found an experiment folder with relative path %s.  Can't synch because that's the path to the to-process folder", ...
                          parent_absolute_path) ;
        else            
            result = {parent_relative_path} ;
        end
    else            
        % For each folder, recurse
        result = cell(0,1) ;
        for i = 1 : length(folder_names) ,
            folder_name = folder_names{i} ;
            child_folder_relative_path_list = ...
                 find_remote_experiment_folders_helper(user_name, ...
                                                       host_name, ...
                                                       fullfile(parent_relative_path, folder_name), ...
                                                       root_absolute_path, ...
                                                       to_process_folder_name, ...
                                                       spinner) ;
            result = [ result ; child_folder_relative_path_list ] ;  %#ok<AGROW>
        end
    end
end
