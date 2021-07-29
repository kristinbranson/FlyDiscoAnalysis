function [relative_path_from_experiment_index, is_aborted_from_experiment_index] = ...
        find_remote_experiment_folders_helper(user_name, host_name, parent_relative_path, root_absolute_path, to_process_folder_name, spinner)
    % Find the experiment folders on a remote host.  Returns relative paths,
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
            is_aborted_from_experiment_index = any(strcmp('ABORTED', file_names)) ;
            relative_path_from_experiment_index = {parent_relative_path} ;
        end
    else            
        % For each folder, recurse
        relative_path_from_experiment_index = cell(0,1) ;
        is_aborted_from_experiment_index = false(0,1) ;
        for i = 1 : length(folder_names) ,
            folder_name = folder_names{i} ;
            [relative_path_from_child_experiment_index, is_aborted_from_child_experiment_index] = ...
                 find_remote_experiment_folders_helper(user_name, ...
                                                       host_name, ...
                                                       fullfile(parent_relative_path, folder_name), ...
                                                       root_absolute_path, ...
                                                       to_process_folder_name, ...
                                                       spinner) ;
            relative_path_from_experiment_index = [ relative_path_from_experiment_index ; relative_path_from_child_experiment_index ] ;  %#ok<AGROW>
            is_aborted_from_experiment_index = [ is_aborted_from_experiment_index ; is_aborted_from_child_experiment_index ] ;  %#ok<AGROW>
        end
    end
end
