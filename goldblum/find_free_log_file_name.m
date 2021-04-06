function result = find_free_log_file_name(experiment_folder_path)
    max_index_to_try = 99 ;
    did_find_free_name = false ;
    for i = 1 : max_index_to_try ,        
        file_name = sprintf('flydisco-analysis-log-run-%02d.txt', i) ;
        result = fullfile(experiment_folder_path, file_name) ;
        if ~exist(result, 'file') ,
            did_find_free_name = true ;
            break
        end
    end
    if ~did_find_free_name ,
        error('Folder %s has too many log files', experiment_folder_path) ;
    end
end
