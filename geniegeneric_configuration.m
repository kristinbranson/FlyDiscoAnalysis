function result = geniegeneric_configuration()
    this_folder_path = fileparts(mfilename('fullpath')) ;
    result = struct() ;
    result.cluster_billing_account_name = 'genie' ;
    result.settings_folder_path = '/groups/genie/home/geniegeneric/GenieFlyDiscoSettings' ;
end
