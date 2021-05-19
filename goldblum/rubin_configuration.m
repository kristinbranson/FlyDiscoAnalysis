function result = rubin_configuration()
    result = struct() ;
    result.lab_head_last_name = 'rubin' ;
    result.host_name_from_rig_index = {'flybowl-ww1.hhmi.org', 'flybowl-ww3.hhmi.org'} ;
    result.rig_user_name_from_rig_index = {'labadmin', 'labadmin'} ;
    result.data_folder_path_from_rig_index = {'/cygdrive/h/flydisco_data', '/cygdrive/e/flydisco_data'} ;
    result.destination_folder = '/groups/rubin/data0/rubinlab/flydisco_data' ;
    this_folder_path = fileparts(mfilename('fullpath')) ;
    flydisco_analysis_path = fileparts(this_folder_path) ;
    result.settings_folder_path = fullfile(flydisco_analysis_path, 'settings') ;
    result.does_have_per_user_folders = true ;
end
