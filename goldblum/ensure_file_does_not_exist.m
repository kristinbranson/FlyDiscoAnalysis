function ensure_file_does_not_exist(raw_file_path)
    file_path = absolute_filename(raw_file_path) ;
    ensure_file_does_not_exist_helper(file_path) ;
end



function ensure_file_does_not_exist_helper(file_path)
    if exist(file_path, 'file') ,
        if exist(file_path, 'dir') ,
            error('Want to delete file at %s, but there''s a folder there, and I''m scared to delete a whole subtree', file_path) ;
        else
            % Presumably a regular file (or a symlink), so we delete it
            delete(file_path)
        end
    else
        % nothing to do, already there's no file at that path
    end        
end
