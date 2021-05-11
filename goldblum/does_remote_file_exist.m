function does_exist = does_remote_file_exist(user_name, host_name, path) 
    escaped_source_path = escape_path_for_bash(path) ;
    remote_stat_command_line = horzcat('test -a ', escaped_source_path) ;
    command_line = horzcat('ssh ', user_name, '@', host_name, ' ', remote_stat_command_line) ;
    [return_code, stdout] = system(command_line) ;
    if return_code==0 ,
        does_exist = true ;
    elseif return_code == 1 ,
        does_exist = false ;
    else
        error('Ambiguous result from "%s": Not clear if file/folder %s exists or not on host %s.  Return code is %d.  stdout is:\n%s', ...
              command_line, path, host_name, return_code, stdout) ;
    end
end
