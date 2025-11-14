function result = realpath(file_name)
    % Get the canonical path (resolving symlinks, '..', etc) from a file name
    if ispc()
      % Half-ass it on Windows
      if is_filename_absolute(file_name)
        result = file_name ;
      else
        result = fullfile(pwd(), file_name) ;
      end
    else
      stdout = system_from_list_with_error_handling({'realpath', file_name}) ;
      result = strtrim(stdout) ;
    end
end
