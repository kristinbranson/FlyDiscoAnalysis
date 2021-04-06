function result = is_an_analysis_protocol_folder(folder_name)
    if exist(folder_name, 'dir') ,
        putative_dataloc_params_txt_file_name = fullfile(folder_name, 'dataloc_params.txt') ;
        result = logical( exist(putative_dataloc_params_txt_file_name, 'file') ) ;
    else
        result = false ;
    end
end
