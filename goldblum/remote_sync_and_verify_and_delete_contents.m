function remote_sync_and_verify_and_delete_contents(source_user, source_host, source_path, dest_path, does_use_per_user_folders)
    remote_sync_and_verify(source_user, source_host, source_path, dest_path) ;
    delete_stuff_in_rig_data_folder(source_user, source_host, source_path, does_use_per_user_folders) ;
end
