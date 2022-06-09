% cd FlyDiscoAnalysis
% matlab &
% <Now in Matlab>
%%
modpath
% lab_head_last_name = 'branson' ;
cluster_billing_account_name = 'branson';
do_use_bqueue = false ;    
do_actually_submit_jobs = false ;  

%% set params                                  
settings_folder_path = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/settings' ;
%%%% this didn't work run with default settings dir based on screen_type
% analysis_protocol = '20210806_flybubble_LED_analysisparams';
% do_force_analysis = true; 
% analysis_parameters = {'analysis_protocol',analysis_protocol};
%%%%
% analysis_protocol = '20211014_flybubbleRed_LED';
% analysis_protocol = '20220217_flybubble_TrpA';
% analysis_protocol = '20220217_flybubble_TrpA_Male';
analysis_protocol = '20220517_flybubble_LED';
analysis_parameters = {'analysis_protocol',analysis_protocol};

% analysis_parameters = {'analysis_protocol',analysis_protocol, ... 
%     'doautomaticchecksincoming','force',...
%     'doflytracking','on', ...
%     'doregistration','on',...
%     'doledonoffdetection','on',...
%     'dosexclassification','on',...
%     'dotrackwings','off',...
%     'docomputeperframefeatures','on',...
%     'docomputehoghofperframefeatures','off',...
%     'dojaabadetect','off',...
%     'docomputeperframestats','off',...
%     'doplotperframestats','off',...
%     'domakectraxresultsmovie','on',...
%     'doextradiagnostics','off',...
%     'doanalysisprotocol',isunix,...
%     'doautomaticcheckscomplete','force'};
% % quick test ACI and ACC
% analysis_parameters = {'analysis_protocol',analysis_protocol, ... 
%     'doautomaticchecksincoming','force',...
%     'doflytracking','off', ...
%     'doregistration','on',...
%     'doledonoffdetection','off',...
%     'dosexclassification','off',...
%     'dotrackwings','off',...
%     'docomputeperframefeatures','off',...
%     'docomputehoghofperframefeatures','off',...
%     'dojaabadetect','off',...
%     'docomputeperframestats','off',...
%     'doplotperframestats','off',...
%     'domakectraxresultsmovie','off',...
%     'doextradiagnostics','off',...
%     'doanalysisprotocol',isunix,...
%     'doautomaticcheckscomplete','force'};

%% make explist
% running as bransonlab 
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210806_testingcaboose/VNC_JRC_SS49220_RigB_20210419T144428', ...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210806_testingcaboose/VNC_JRC_SS49220_RigB_20210421T143507'};

% make explist of only experiment dirs with tracking and NOT aborted
% rootdatadir = '/groups/branson/bransonlab/flydisco_data'; 
% lab ='branson';
% outdir4list = '/groups/branson/home/robiea/Projects_data/FlyDisco/FlyDiscoPipeline';

% rootdatadir = '/groups/dickson/dicksonlab/flydisco_data';
% lab = 'dickson';
% outdir4list ='/groups/dickson/dicksonlab/Alice';
% 
% filesavestr = ['RERUNposttracking_explist_',lab,datestr(now,'yyyymmddTHHMMSS'),'.txt'];
% filesavename = fullfile(outdir4list,filesavestr);
% fid = fopen(filesavename,'w');
% 
% trackerfilestr = 'movie-track.mat';
% metadatafile = 'Metadata.xml';
% screen_type = 'non_olympiad_dickson_VNC';
% 
% % pull all the data with VNC* screen_type 
% % bransonlab
% 
% 
% 
% [expdirstruct] = getExperimentDirsFlyDisco(rootdatadir,'screen_type','VNC*');
% metadata = expdirstruct;
% 
% % changed to NO tracking files for running expdirs that didn't process at
% % all
% for i = 1:numel(metadata)
%     expdir = metadata(i).file_system_path;    
%     trackerfile = fullfile(expdir,trackerfilestr);
%     abortedfile = fullfile(expdir,'ABORTED');    
%     if ~exist(trackerfile,'file') && ~exist(abortedfile,'file') && strcmp(metadata(i).screen_type,screen_type)
%         fprintf(fid,'%s\n',expdir);
%     end
% end
% 
% 
% folder_path_from_experiment_index = textread(filesavename,'%s');

%for testing on first 5 experiments

% folder_path_from_experiment_index = folder_path_from_experiment_index(1);
%% load explist for rerunning caboose
% folder_path_from_experiment_index = textread('/groups/dickson/dicksonlab/Alice/RERUNposttracking_explist_dickson20210809T141224.txt','%s');
% folder_path_from_experiment_index = folder_path_from_experiment_index(1);

% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20210913_testingNorpA/CsChr_JHS_K_85321_RigA_20210903T085205'};

% run pipeline for jaaba file conversion
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/JAABA/Data_FlyBubble/FlyTracker/cx_JHS_K_85321_CsChr_RigD_20150909T163219_tmp'};
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/JAABA/Data_FlyBubble/FlyTracker/cx_GMR_SS00030_CsChr_RigC_20150826T144616',...
% '/groups/branson/home/robiea/Projects_data/JAABA/Data_FlyBubble/FlyTracker/cx_GMR_SS00038_CsChr_RigB_20150729T150617'};
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/JAABA/Data_FlyBubble/FlyTracker/cx_GMR_SS00168_CsChr_RigD_20150909T111218'};

% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/test_20220217_flybubble_TrpA/TrpAFemale2_GMR_72C11_vk5_RigA_20220216T090844'};
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/test_20220217_flybubble_TrpA/TrpAFemale2_GMR_Eb5_vk5_RigB_20220216T085536',...
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/test_20220217_flybubble_TrpA/TrpAMale2_GMR_71G01_JK73A_RigB_20220216T074837'};
%test movie params for 20220414 LED protocol
% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/test_20220517_flybubble_VNC2_testmovieparams/VNC2_EXT_VGLUT-GAL4_RigA_20220511T094028',...
% '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/test_20220517_flybubble_VNC2_testmovieparams/VNC2_YNA_K_162984_RigD_20220511T103822'};
% test rerun 6/9/2022 gender wrong
folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/test_20220608_testpipelinefiledeletion/VNC2_JRC_SS83407_RigB_20220419T092554'};

%% delete pipeline files before rerunning pipeline
% use CleanOutExpDirs or 
CleanOutExpDirs_leavetracking(folder_path_from_experiment_index)


%% run analysis
% SELECTED analysis_parameters
% goldblum_analyze_experiment_folders(folder_path_from_experiment_index, settings_folder_path, cluster_billing_account_name, ...
%     do_use_bqueue, do_actually_submit_jobs, analysis_parameters) ;

% DEFAULT analysis_parameters
goldblum_analyze_experiment_folders(folder_path_from_experiment_index, settings_folder_path, cluster_billing_account_name, ...
    do_use_bqueue, do_actually_submit_jobs) ;