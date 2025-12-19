function result = reiserrobot_configuration()
    result = struct() ;
    result.cluster_billing_account_name = 'reiser' ;
    % fda_folder_path = fileparts(fileparts(mfilename('fullpath'))) ;
    % result.settings_folder_path = fullfile(fda_folder_path, 'settings-internal') ;
    result.settings_folder_path = '/groups/reiser/home/reiserrobot/ReiserFlyDiscoSettings/settings' ;
end
