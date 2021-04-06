function delete_remote_folder_contents(remote_user_name, remote_dns_name, remote_folder_path) 
    % Delete the contents of the remote folder after successful transfer
    escaped_remote_folder_path = escape_path_for_bash(remote_folder_path) ;
    remote_command_line = ...
        sprintf('ls -1 -A %s | while read line; do echo %s/$line; done | xargs -d ''\\n'' rm -rf', escaped_remote_folder_path, escaped_remote_folder_path) ;
    escaped_remote_command_line = escape_path_for_bash(remote_command_line) ;
    command_line = sprintf('ssh %s@%s %s', remote_user_name, remote_dns_name, escaped_remote_command_line) ;
    [return_code, stdout] = system(command_line) ;
    if return_code ~= 0 ,
        error('Unable to delete the contents of folder %s as user %s on host %s.  Stdout/stderr was:\n%s\n', ...
                remote_folder_path, remote_user_name, remote_dns_name, stdout) ;
    end
end
