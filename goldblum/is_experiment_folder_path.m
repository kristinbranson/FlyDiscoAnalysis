function result = is_experiment_folder_path(path) 
    putative_ufmf_file_path = fullfile(path, 'movie.ufmf') ;
    if exist(putative_ufmf_file_path, 'file') ,
        result = true ;
        return
    end
    putative_avi_file_path = fullfile(path, 'movie.avi') ;
    if exist(putative_avi_file_path, 'file') ,
        result = true ;
        return
    end    
    result = false ;
end
