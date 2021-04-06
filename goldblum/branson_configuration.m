function result = branson_configuration()
    result = struct() ;
    result.lab_head_last_name = 'branson' ;
    result.host_name_from_rig_index = { 'arrowroot.hhmi.org', 'beet.hhmi.org', 'carrot.hhmi.org', 'daikon.hhmi.org' } ;
    result.rig_user_name_from_rig_index = repmat({'bransonk'}, [1 4]) ;
    result.data_folder_path_from_rig_index = repmat({'/cygdrive/e/flydisco_data'}, [1 4]) ;
    result.destination_folder = '/groups/branson/bransonlab/flydisco_data' ;
    this_folder_path = fileparts(mfilename('fullpath')) ;
    result.settings_folder_path = fullfile(this_folder_path, 'FlyDiscoAnalysis/settings') ;
    result.does_use_per_user_folders = false ;
end
