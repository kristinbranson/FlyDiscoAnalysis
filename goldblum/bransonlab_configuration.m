function result = bransonlab_configuration()
    result = struct() ;
    result.cluster_billing_account_name = 'branson' ;
    result.host_name_from_rig_index = { 'arrowroot.hhmi.org', 'beet.hhmi.org', 'carrot.hhmi.org', 'daikon.hhmi.org' } ;
    result.rig_user_name_from_rig_index = repmat({'bransonk'}, [1 4]) ;
    result.data_folder_path_from_rig_index = repmat({'/cygdrive/e/flydisco_data/branson'}, [1 4]) ;
    result.destination_folder = '/groups/branson/bransonlab/flydisco_data' ;
    %this_folder_path = fileparts(mfilename('fullpath')) ;
    %flydisco_analysis_path = fileparts(this_folder_path) ;
    %result.settings_folder_path = fullfile(flydisco_analysis_path, 'settings') ;
    result.settings_folder_path = '/groups/branson/home/bransonlab/settings' ;
    result.does_use_per_user_folders = false ;
end
