function delete_remote_folder(remote_user_name, remote_dns_name, remote_folder_path)
    escaped_remote_folder_path = escape_for_bash(remote_folder_path) ;
    command_line = sprintf('ssh %s@%s rm -rf %s', remote_user_name, remote_dns_name, escaped_remote_folder_path) ;
    return_code = system(command_line) ;
    if return_code ~= 0 ,
        error('Unable to delete the folder %s as user %s on host %s', ...
              remote_folder_path, remote_user_name, remote_dns_name) ;
    end
end
    