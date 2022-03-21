function push_goldblum_into_production_for_branson_lab()
    % Determine the FlyDiscoAnalysis folder path
    goldblum_folder_path = fileparts(mfilename('fullpath')) ;
    fda_folder_path = fileparts(goldblum_folder_path) ;
    
    % Make sure there are no uncommitted changes
    error_if_uncommited_changes(fda_folder_path) ;
    
    % Do Branson Lab instance
    lab_name = 'bransonlab' ;
    copy_to_single_user_account(lab_name, fda_folder_path) ;    

    % If get here, everything went well
    fprintf('Successfully copied %s into the %s user account\n', fda_folder_path, lab_name) ;
end



function copy_to_single_user_account(user_name, fda_folder_path)
    % Copy the folder over
    host_name = 'login2' ;  % Why not?
    fprintf('Copying into the %s user account...', user_name) ;
    remote_system_from_list_with_error_handling(user_name, host_name, {'rm', '-rf', 'FlyDiscoAnalysis'}) ;
    remote_system_from_list_with_error_handling(user_name, host_name, {'cp', '-R', '-T', fda_folder_path, 'FlyDiscoAnalysis'}) ;
    fprintf('done.\n') ;
end
