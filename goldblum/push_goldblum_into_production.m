function push_goldblum_into_production(account_name_from_lab_index)
    % Deal with arguments
    if ~exist('account_name_from_lab_index', 'var') || isempty(account_name_from_lab_index) ,
        % If no arg, do all labs
        account_name_from_lab_index = { 'bransonlab', 'rubinlab', 'projtechreslab', 'geniegeneric' } ;
    end
    
    % If char arg, convert to cellstring
    if ischar(account_name_from_lab_index) ,
        account_name_from_lab_index = { account_name_from_lab_index } ;
    end
    
    % Determine the FlyDiscoAnalysis folder path
    goldblum_folder_path = fileparts(mfilename('fullpath')) ;
    fda_folder_path = fileparts(goldblum_folder_path) ;
    
    % Make sure there are no uncommitted changes
    error_if_uncommited_changes(fda_folder_path) ;
    
    % Do each lab
    cellfun(@(account_name)(copy_to_single_user_account(account_name, fda_folder_path)), ...
            account_name_from_lab_index) ;

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
