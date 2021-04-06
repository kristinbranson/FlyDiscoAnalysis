function result = find_experiment_folders_helper(parent_path, spinner)
    % Get a list of all files and folders in the source, dest folders
    [entries, is_entry_a_folder] = simple_dir(parent_path) ;
    spinner.spin() ;

    % Separate source file, folder names    
    file_names = entries(~is_entry_a_folder) ;
    folder_names = entries(is_entry_a_folder) ;    

    % If the parent_path is an experiment folder, we're done
    if is_experiment_folder_given_contents(file_names) ,
        result = {parent_path} ;
    else            
        % For each folder, recurse
        result = cell(0,1) ;
        for i = 1 : length(folder_names) ,
            folder_name = folder_names{i} ;
            child_folder_path_list = ...
                 find_experiment_folders_helper(fullfile(parent_path, folder_name), ...
                                                spinner) ;
            result = [ result ; child_folder_path_list ] ;  %#ok<AGROW>
        end
    end
end
