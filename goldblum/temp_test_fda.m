experiment_folder_path = '/groups/branson/bransonlab/taylora/flydisco/FlyDiscoAnalysis/goldblum/goldblum-test-destination-folder/SS36564_20XUAS_CsChrimson_mVenus_attP18_flyBowlMing_20200227_Continuous_2min_5int_20200107_20200229T132141'
settings_folder_path = '/groups/branson/bransonlab/taylora/flydisco/FlyDiscoAnalysis/settings' 

FlyDiscoPipeline(experiment_folder_path, struct('settingsdir', {settings_folder_path})) ;
