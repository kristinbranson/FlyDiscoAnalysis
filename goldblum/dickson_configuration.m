function result = dickson_configuration()
    result = struct() ;
    result.lab_head_last_name = 'dickson' ;
    result.host_name_from_rig_index = { 'arrowroot.hhmi.org', 'beet.hhmi.org', 'carrot.hhmi.org', 'daikon.hhmi.org' } ;
    result.rig_user_name_from_rig_index = repmat({'bransonk'}, [1 4]) ;
    result.data_folder_path_from_rig_index = repmat({'/cygdrive/e/flydisco_data'}, [1 4]) ;
    result.destination_folder = '/groups/dickson/dicksonlab/flydisco_data' ;
    this_folder_path = fileparts(mfilename('fullpath')) ;
    flydisco_analysis_path = fileparts(this_folder_path) ;
    result.settings_folder_path = fullfile(flydisco_analysis_path, 'settings') ;
    result.does_use_per_user_folders = false ;
end
