function result = is_experiment_folder_given_name(folder_name)
    % Given the name (or path) to a folder, determine whether the folder is a
    % FlyDisco experiment folder.  Returns a logical scalar that is true iff the
    % folder is an experiment folder.

    % Get the file names in the folder
    [entry_names, is_entry_a_folder] = simple_dir(folder_name) ;
    file_names = entry_names(~is_entry_a_folder) ;

    % Call the function that makes the determination, based on file names
    result = is_experiment_folder_given_contents(file_names) ;
end
