
expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240508_VNC3/VNC3_YNA_K_162984_RigA_20240507T182023';
settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis';
analysis_protocol = 'settings-internal/20240507_flybubble_LED_VNC3';
load('/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20240508_VNC3/VNC3_YNA_K_162984_RigA_20240507T182023/BAK_indicatordata.mat','indicatorLED');

error_if_protocol_stim_num_notequal_detected(expdir, settingsdir, analysis_protocol,'indicatorLED',indicatorLED)