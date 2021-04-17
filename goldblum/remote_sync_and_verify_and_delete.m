function remote_sync_and_verify_and_delete(source_user_name, source_host_name, source_path, dest_path) 
    % Synch the destination folder to the source folder
    remote_sync(source_user_name, source_host_name, source_path, dest_path) ;

    % Now verify that that worked (will raise an exception if verification fails)
    remote_verify(source_user_name, source_host_name, source_path, dest_path) ;
    
    % Finally, delete the remote folder
    delete_remote_folder(source_user_name, source_host_name, source_path) ;
end
