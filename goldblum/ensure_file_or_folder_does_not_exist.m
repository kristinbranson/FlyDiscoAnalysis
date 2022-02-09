function ensure_file_or_folder_does_not_exist(raw_file_path)
    file_path = absolute_filename(raw_file_path) ;
    ensure_file_or_folder_does_not_exist_helper(file_path) ;
end



function ensure_file_or_folder_does_not_exist_helper(file_path)
    if exist(file_path, 'file') ,
        if exist(file_path, 'dir') ,
            system_from_list_with_error_handling({'rm', '-rf', file_path}) ;
        else
            % Presumably a regular file (or a symlink), so we delete it
            delete(file_path)
        end
    else
        % nothing to do, already there's no file at that path
    end        
end
