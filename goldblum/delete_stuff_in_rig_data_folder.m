function delete_stuff_in_rig_data_folder(remote_user_name, remote_dns_name, root_remote_folder_path, does_use_per_user_folders)
    % Delete the contents of the remote folder after successful transfer
    if does_use_per_user_folders ,
        try
            [remote_entries, ~, ~, is_remote_entry_a_dir] = ...
                list_remote_dir(remote_user_name, remote_dns_name, root_remote_folder_path) ;
        catch me ,
            % orginally caught OSError
            if isequal(me.identifier, 'list_remote_dir:failed') ,
                fprintf('Unable to delete the contents of folder %s as user %s on host %s\n', root_remote_folder_path, remote_user_name, remote_dns_name) ;
                disp(me.getReport()) ;
            else
                rethrow(me) ;
            end
        end

        % separate remote dirs out (these should be the user folders)
        remote_user_folder_names = remote_entries(is_remote_entry_a_dir) ;

        % delete the contents of each of the remote user folders
        for i = 1 : length(remote_user_folder_names) ,
            remote_user_folder_name = remote_user_folder_names{i} ;
            remote_user_folder_path = fullfile(root_remote_folder_path, remote_user_folder_name) ;
            delete_remote_folder_contents(remote_user_name, remote_dns_name, remote_user_folder_path) ;
        end
    else
        % If not using per-user folders, can just delete contents of whole
        % folder
        delete_remote_folder_contents(remote_user_name, remote_dns_name, root_remote_folder_path) ;
    end
end
