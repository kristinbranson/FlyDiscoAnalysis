function result = otopalikrobot_configuration()
    this_folder_path = fileparts(mfilename('fullpath')) ;
    result = struct() ;
    result.cluster_billing_account_name = 'otopalik' ;
    result.settings_folder_path = fullfile(this_folder_path, 'settings-internal') ;
end
