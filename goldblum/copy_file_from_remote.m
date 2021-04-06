function elapsed_time = copy_file_from_remote(source_user, source_host, source_path, dest_path)
    escaped_source_path = escape_path_for_bash(source_path) ;
    escaped_dest_path = escape_path_for_bash(dest_path) ;
    source_spec = horzcat(source_user, '@', source_host, ':', escaped_source_path) ;
    tic_id = tic() ;
    command_line = sprintf('scp -v -B -T %s %s', source_spec, escaped_dest_path) ;
    [scp_return_code, stdout] = system(command_line) ;
    % scp doesn't honor the user's umask, so we need to set the file
    % permissions explicitly
    if scp_return_code == 0,
        command_line = sprintf('chmod u+rw-x,g+rw-x,o+r-wx %s', escaped_dest_path) ;
        [chmod_return_code, stdout] = system(command_line) ;
    end
    elapsed_time = toc(tic_id) ;
    if scp_return_code ~= 0 ,
        error('copy_file_from_remote:failed', ...
              'Unable to copy the file %s as remote user %s from host %s to destination %s:\n%s', ... 
              source_path, source_user, source_host, dest_path, stdout) ;
    elseif chmod_return_code ~= 0 ,
        error('copy_file_from_remote:failed', ...
              'Unable to set the permissions of %s after copy:\n%s', ...
              dest_path, stdout) ;
    end    
end
