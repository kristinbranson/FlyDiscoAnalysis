function metadata_file_path = determine_metadata_file_path(experiment_folder_path)
    [name_from_entry_index, is_folder_from_entry_index] = simple_dir(experiment_folder_path) ;
    is_file_from_entry_index = ~is_folder_from_entry_index ;
    name_from_file_index = name_from_entry_index(is_file_from_entry_index) ;
    is_metadata_file = strcmpi('metadata.xml', name_from_file_index) ;
    name_from_metadata_file_index = name_from_file_index(is_metadata_file) ;
    if isempty(name_from_metadata_file_index) ,
        error('There seem to be zero files in %s having names that match "metadata.xml" when ignoring case', experiment_folder_path) ;
    elseif isscalar(name_from_metadata_file_index) ,
        metadata_file_name = name_from_metadata_file_index{1} ;
        metadata_file_path = fullfile(experiment_folder_path, metadata_file_name) ;
    else
        error('There seem to be %d files in %s having names that match "metadata.xml" when ignoring case', length(name_from_metadata_file_index)) ;
    end
end
