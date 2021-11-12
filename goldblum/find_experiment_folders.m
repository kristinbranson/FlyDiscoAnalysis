function result = find_experiment_folders(source_path)
    relative_path_from_experiment_index = find_experiment_folders_relative(source_path) ;
    result = cellfun(@(rel_path)(fullfile(source_path, rel_path)), relative_path_from_experiment_index, 'UniformOutput', false) ;
end
