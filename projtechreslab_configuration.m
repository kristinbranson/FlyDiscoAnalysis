function result = projtechreslab_configuration()
    result = struct() ;
    result.cluster_billing_account_name = 'projtechres' ;
    result.host_name_from_rig_index = {'flybowl-ww1.hhmi.org'} ;
    result.rig_user_name_from_rig_index = {'labadmin'} ;
    result.data_folder_path_from_rig_index = {'/cygdrive/e/flydisco_data/projtechres'} ;
    result.destination_folder = '/groups/projtechres/projtechres/flydisco_data' ;
    this_folder_path = fileparts(mfilename('fullpath')) ;
    result.settings_folder_path = fullfile(this_folder_path, 'settings-internal') ;
    result.does_have_per_user_folders = true ;
end
