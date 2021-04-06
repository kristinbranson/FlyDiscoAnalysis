function [experiment_dir_names, other_dir_names] = separate_experiment_folders_from_others(parent_path, source_dir_names) 
    is_experiment_dir = false(size(source_dir_names)) ;
    for i = 1 : length(source_dir_names) ,
        source_dir_name = source_dir_names{i} ;
        is_experiment_dir(i) = is_experiment_folder_path(fullfile(parent_path, source_dir_name)) ;
    end
    experiment_dir_names = source_dir_names(is_experiment_dir) ;
    other_dir_names = source_dir_names(~is_experiment_dir) ;
end
