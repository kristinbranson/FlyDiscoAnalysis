% cd FlyDiscoAnalysis
% matlab &
% <Now in Matlab>
%%
modpath
lab_head_last_name = 'branson' ;
do_use_bqueue = false ;    
do_actually_submit_jobs = false ;  

%% set params                                  
settings_folder_path = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/settings' ;

% analysis_protocol = '20210319_flybubble_LED_testing';
% analysis_protocol = 'current_non_olympiad_dickson_VNC';
% folder_path_from_experiment_index = { '/groups/branson/bransonlab/flydisco_data/VNC_YNA_K_162984_RigD_20210405T152704' } ;
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210407_testingpipeline/VNC_YNA_K_162984_RigD_20210405T152704'};
% test with Kristin's Nan removal
% analysis_protocol = '20210409_flybubble_nonLED';
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210409_testingpipeline/nochr_TrpA71G01_Unknown_RigA_20201216T162938'};
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210409_testingpipeline/nochr_TrpA65F12_Unknown_RigD_20201216T175902'};
% analysis_protocol = '20210319_flybubble_LED_testing';

% analysis_protocol = 'current_non_olympiad_dickson_VNC';
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/TestSuite/FTtracked/FlyBubbleRGB/LED/VNC_JHS_K_85321_RigA_20210408T130721'};
analysis_protocol =  '20210417_flybubble_LED';
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/TestSuite/FTtracked/FlyBubbleRGB/LED/VNC_JHS_K_85321_RigB_20210408T130659'};
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/TestSuite/FTtracked/FlyBubbleRGB/LED/VNC_JHS_K_85321_RigC_20210408T124720'};
folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/TestSuite/FTtracked/FlyBubbleRGB/LED/VNC_JHS_K_85321_RigD_20210408T130906'};
% whether or not to use fial/sucess flags
do_force_analysis = true;    

% analysis_parameters = cell(1, 0) ;

analysis_parameters = {'analysis_protocol',analysis_protocol, ...    
    'doautomaticchecksincoming','on',...
    'doflytracking','on', ...
    'doregistration','on',...
    'doledonoffdetection','on',...
    'dosexclassification','on',...
    'dotrackwings','off',...
    'docomputeperframefeatures','on',...
    'docomputehoghofperframefeatures','off',...
    'dojaabadetect','off',...
    'docomputeperframestats','off',...
    'doplotperframestats','off',...
    'domakectraxresultsmovie','on',...
    'doextradiagnostics','off',...
    'doanalysisprotocol',isunix,...
    'doautomaticcheckscomplete','off'};

%% run analysis

analyze_experiment_folders(folder_path_from_experiment_index, settings_folder_path, lab_head_last_name, ...
                           do_force_analysis, do_use_bqueue, do_actually_submit_jobs, analysis_parameters) ;