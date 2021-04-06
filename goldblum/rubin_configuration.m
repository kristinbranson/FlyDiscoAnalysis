function result = rubin_configuration()
    result = struct() ;
    result.lab_head_last_name = 'rubin' ;
    result.host_name_from_rig_index = {'flybowl-ww1.hhmi.org'} ;
    result.rig_user_name_from_rig_index = {'labadmin'} ;
    result.data_folder_path_from_rig_index = {'/cygdrive/h/flydisco_data'} ;
    result.destination_folder = '/groups/rubin/data0/rubinlab/flybowlanalysis-flytracker-data' ;    
    result.settings_folder_path = '/groups/rubin/home/rubinlabuser/flybowlanalysis-flytracker-settings' ;
    result.does_have_per_user_folders = true ;
end
