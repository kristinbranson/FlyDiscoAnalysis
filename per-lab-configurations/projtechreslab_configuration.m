function result = projtechreslab_configuration()
    fda_folder_path = fileparts(fileparts(mfilename('fullpath'))) ;
    result = struct() ;
    result.cluster_billing_account_name = 'projtechres' ;
    result.settings_folder_path = fullfile(fda_folder_path, 'settings-internal') ;
end
