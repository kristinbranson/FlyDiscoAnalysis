function result = realpath(file_name)
    % Get the canonical path (resolving symlinks, '..', etc) from a file name    
    stdout = system_from_list_with_error_handling({'realpath', file_name}) ;
    result = strtrim(stdout) ;
end
