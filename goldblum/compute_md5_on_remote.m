function hex_digest =  compute_md5_on_remote(source_user, source_host, source_path)
    escaped_source_path = escape_path_for_bash(source_path) ;
    host_spec = horzcat(source_user, '@', source_host) ;
    command_line = sprintf('ssh %s md5sum %s', host_spec, escaped_source_path) ;
    [return_code, stdout] = system(command_line) ;
    if return_code ~= 0 ,
        error('Unable to md5sum the file %s as user %s on host %s', source_path, source_user, source_host) ;
    end
    tokens = strsplit(stdout) ;
    if isempty(tokens) ,
        error('Got a weird result while md5sum''ing the file %s as user %s on host %s.  Result was: %s', ...
              source_path, source_user, source_host, stdout) ;
    end
    hex_digest = tokens{1} ;
end
    