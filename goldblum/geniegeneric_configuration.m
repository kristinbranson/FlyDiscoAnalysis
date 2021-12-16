function result = geniegeneric_configuration()
    result = struct() ;
    result.cluster_billing_account_name = 'genie' ;
    result.host_name_from_rig_index = {'flybowl-ww1.hhmi.org'} ;
    result.rig_user_name_from_rig_index = {'labadmin'} ;
    result.data_folder_path_from_rig_index = {'/cygdrive/e/flydisco_data/genie'} ;
    result.destination_folder = '/groups/genie/genie/flydisco_data' ;
    this_folder_path = fileparts(mfilename('fullpath')) ;
    flydisco_analysis_path = fileparts(this_folder_path) ;
    result.settings_folder_path = fullfile(flydisco_analysis_path, 'settings') ;
    result.does_have_per_user_folders = true ;
end
