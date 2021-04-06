function remote_sync_and_verify(source_user, source_host, source_path, dest_path) 
    % Call the function that does the real work
    remote_sync(source_user, source_host, source_path, dest_path) ;

    % Now verify that that worked (will raise an exception if verification fails)
    remote_verify(source_user, source_host, source_path, dest_path) ;
end
