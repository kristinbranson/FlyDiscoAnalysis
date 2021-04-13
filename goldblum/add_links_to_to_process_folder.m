function add_links_to_to_process_folder(destination_folder, to_process_folder_name, relative_path_from_experiment_folder_index)
    to_process_folder_path = fullfile(destination_folder, to_process_folder_name) ;
    escaped_to_process_folder_path = escape_string_for_bash(to_process_folder_path) ;
    experiment_folder_relative_path_count = length(relative_path_from_experiment_folder_index) ;
    for i = 1 : experiment_folder_relative_path_count ,
        experiment_folder_relative_path = relative_path_from_experiment_folder_index{i} ;
        experiment_folder_absolute_path = fullfile(destination_folder, experiment_folder_relative_path) ;
        escaped_experiment_folder_absolute_path = escape_string_for_bash(experiment_folder_absolute_path) ;
        command_line = sprintf('ln -s %s %s', escaped_experiment_folder_absolute_path, escaped_to_process_folder_path) ;
        system_with_error_handling(command_line) ;
    end
end
