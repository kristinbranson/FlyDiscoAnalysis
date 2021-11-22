function result = find_experiment_folders_relative_helper(root_path, parent_relative_path, spinner)
    % Get a list of all files and folders in the source, dest folders
    parent_path = fullfile(root_path, parent_relative_path) ;
    [entries, is_entry_a_folder] = simple_dir(parent_path) ;
    spinner.spin() ;

    % Separate source file, folder names    
    file_names = entries(~is_entry_a_folder) ;
    raw_folder_names = entries(is_entry_a_folder) ;    
    
    % Exclude the to-process folder if in the root path
    if isempty(parent_relative_path) ,
        is_to_process_folder = strcmp(raw_folder_names, 'to-process') ;
        folder_names = raw_folder_names(~is_to_process_folder) ;
    else
        folder_names = raw_folder_names ;
    end
    
    % If the parent_path is an experiment folder, we're done
    if is_experiment_folder_given_contents(file_names) ,
        result = {parent_relative_path} ;
    else            
        % For each folder, recurse
        result = cell(0,1) ;
        for i = 1 : length(folder_names) ,
            folder_name = folder_names{i} ;
            child_folder_path_list = ...
                 find_experiment_folders_relative_helper(root_path, ...
                                                         fullfile(parent_relative_path, folder_name), ...
                                                         spinner) ;
            result = [ result ; child_folder_path_list ] ;  %#ok<AGROW>
        end
    end
end
