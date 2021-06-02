function copy_goldblum_into_production()
    % Determine the FlyDiscoAnalysis folder path
    goldblum_folder_path = fileparts(mfilename('fullpath')) ;
    fda_folder_path = fileparts(goldblum_folder_path) ;
    
    % Make sure there are no uncommitted changes
    error_if_uncommited_changes(fda_folder_path) ;
    
    % Do Branson Lab instance
    copy_to_single_user_account('bransonlab', fda_folder_path) ;
    
    % Do Dickson Lab instance
    copy_to_single_user_account('dicksonlab', fda_folder_path) ;
    
    % Do Rubin Lab instance
    copy_to_single_user_account('rubinlab', fda_folder_path) ;

    % If get here, everything went well
    fprintf('Successfully copied %s into all the *lab user accounts\n', fda_folder_path) ;
end



function copy_to_single_user_account(user_name, fda_folder_path)
    % Copy the folder over
    host_name = 'login2' ;  % Why not?
    fprintf('Copying into the %s user account...', user_name) ;
    remote_system_from_list_with_error_handling(user_name, host_name, {'rm', '-rf', 'FlyDiscoAnalysis'}) ;
    remote_system_from_list_with_error_handling(user_name, host_name, {'cp', '-R', '-T', fda_folder_path, 'FlyDiscoAnalysis'}) ;
    fprintf('done.\n') ;
end
