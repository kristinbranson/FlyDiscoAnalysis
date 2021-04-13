function delete_remote_folder(remote_user_name, remote_dns_name, remote_folder_path)
    % Delete (as in rm -rf) a remote folder
    
    % Check for a truly horrible mistake
    trimmed_remote_folder_path = strtrim(remote_folder_path) ;
    if isempty(trimmed_remote_folder_path) || isequal(trimmed_remote_folder_path, '/') ,
        error('Yeah, I''m not going to rm -rf / on %s', remote_dns_name) ;
    end
    
    escaped_remote_folder_path = escape_string_for_bash(remote_folder_path) ;
    command_line = sprintf('ssh %s@%s rm -rf %s', remote_user_name, remote_dns_name, escaped_remote_folder_path) ;
    return_code = system(command_line) ;
    if return_code ~= 0 ,
        error('Unable to delete the folder %s as user %s on host %s', ...
              remote_folder_path, remote_user_name, remote_dns_name) ;
    end
end
    