function hex_digest = compute_md5_on_local(local_path) 
    escaped_local_path = escape_path_for_bash(local_path) ;
    command_line = sprintf('md5sum %s', escaped_local_path) ;
    [return_code, stdout] = system(command_line) ;
    if return_code ~= 0 ,
        error('Unable to md5sum the file %s.  Stdout/stderr was:\n%s', local_path, stdout) ;
    end
    tokens = strsplit(strtrim(stdout)) ;
    if isempty(tokens) , 
        error('Got a weird result while md5sum''ing the file %s.  Stdout/stderr was:\n%s', local_path, stdout)
    end        
    hex_digest = tokens{1} ;
end
