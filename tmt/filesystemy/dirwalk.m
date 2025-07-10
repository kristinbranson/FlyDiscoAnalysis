function result = dirwalk(root_folder_name, callback, initial_state, varargin)
    % Like Python's os.walk()
    %
    % callback should be a function that can be called like so:
    % updated_state = callback(root_folder_absolute_path, current_folder_relative_path, folder_names, file_names, state, varargin)
    % where root_folder_absolute_path is the absolute path to root_folder_name,     
    % current_folder_relative_path is a folder under root_folder_name, and
    % folder_names are the folders within it, and file_names are the files within
    % it.  The callback will be called for root_folder_name and for all folders
    % under it, in some unspecified order.  The state input and updated_state output
    % allow information to be accumulated across the callback calls.
    
    root_folder_path = absolute_filename(root_folder_name) ;
    result = dirwalk_helper(root_folder_path, '', callback, initial_state, varargin{:}) ;
end



function state = dirwalk_helper(root_folder_path, current_folder_relative_path, callback, initial_state, varargin)
    current_folder_absolute_path = fullfile(root_folder_path, current_folder_relative_path) ;
    [entries, is_entry_a_folder] = simple_dir(current_folder_absolute_path) ;
    folder_names = entries(is_entry_a_folder) ;
    file_names = entries(~is_entry_a_folder) ;
    state = feval(callback, root_folder_path, current_folder_relative_path, folder_names, file_names, initial_state, varargin{:}) ;
    for i = 1 : length(folder_names) ,
        folder_name = folder_names{i} ;
        state = dirwalk_helper(root_folder_path, fullfile(current_folder_relative_path, folder_name), callback, state, varargin{:}) ;
    end
end
