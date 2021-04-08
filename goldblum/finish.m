function finish()
    try
        % Clean up any host-and-PID-specific scratch folder
        host_name = get_short_host_name() ;
        parpool_data_location_folder_name = ['host-' host_name '-pid-' num2str(feature('getpid'))] ;
        parpool_data_location_folder_path = fullfile(get_scratch_folder_path(), parpool_data_location_folder_name) ;
        if exist(parpool_data_location_folder_path, 'file') ,
            escaped_parpool_data_location_folder_path = escape_string_for_bash(parpool_data_location_folder_path) ;            
            command_line = sprintf('rm -rf %s', escaped_parpool_data_location_folder_path) ;
            system_with_error_handling(command_line) ;
        end
    catch me
        fprintf('There was a problem during execution of the finish() function:\n') ;
        fprintf('%s\n', me.getReport()) ;
    end
end
