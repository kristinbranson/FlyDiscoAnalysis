function result = scicompsoft_configuration()
    result = struct() ;
    result.cluster_billing_account_name = 'scicompsoft' ;
    result.host_name_from_rig_index = {'flybowl-ww1.hhmi.org'} ;
    result.rig_user_name_from_rig_index = {'labadmin'} ;
    result.data_folder_path_from_rig_index = {'/cygdrive/h/flydisco_data/scicompsoft'} ;
    result.destination_folder = '/groups/branson/bransonlab/taylora/flydisco/goldblum/goldblum-test-destination-folder' ;    
    result.settings_folder_path = '/groups/branson/bransonlab/taylora/flydisco/goldblum/FlyDiscoAnalysis/settings' ;
    result.does_use_per_user_folders = true ;
end
