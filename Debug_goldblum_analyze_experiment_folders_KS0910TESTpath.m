% cd FlyDiscoAnalysis
% matlab &
% <Now in Matlab>
%%
modpath
lab_head_last_name = 'rubin' ;
do_use_bqueue = false ;    
do_actually_submit_jobs = false ;  

%% set params                                  
%settings_folder_path = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/settings' ;
settings_folder_path = '/groups/rubin/home/schretterc/Documents/FlyDiscoAnalysis/FlyDiscoAnalysis-new/settings' ;
%%%% this didn't work run with default settings dir based on screen_type
% analysis_protocol = '20210806_flybubble_LED_analysisparams';
% do_force_analysis = true; 
%analysis_parameters = {'analysis_protocol',analysis_protocol};
%%%%
%analysis_protocol = '20210909_FlyBowlRGBVision_addperframe';
do_force_analysis = true;   
analysis_parameters = {'doautomaticchecksincoming','on', ...
    'doflytracking','on', ...
    'doregistration','on',...
    'doledonoffdetection','on',...
    'dosexclassification','on',...
    'dotrackwings','on',...
    'docomputeperframefeatures','on',...
    'docomputehoghofperframefeatures','off',...
    'dojaabadetect','off',...
    'docomputeperframestats','off',...
    'doplotperframestats','off',...
    'domakectraxresultsmovie','on',...
    'doextradiagnostics','off',...
    'doanalysisprotocol',isunix,...
    'doautomaticcheckscomplete','on'};
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
folder_path_from_experiment_index = {'/groups/rubin/home/schretterc/Documents/FlyDiscoAnalysis_ExptsToAnalyze_Test/20210401T134552_rig1_flyBowl4__aIPgSS1UASCsChrimson_KS_redonly_protocolRGB_0315_2'};
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
%folder_path_from_experiment_index = textread('/groups/dickson/dicksonlab/Alice/RERUNposttracking_explist_dickson20210809T141224.txt','%s');
%experiment_folder_path = '/groups/rubin/home/schretterc/Documents/FlyDiscoAnalysis_ExptsToAnalyze_Test/20210819T085737_rig1_flyBowl3__20XUASCsChrimsonattp18_SS36564_KS_redonly3times10AND30_080521';
folder_path_from_experiment_index = folder_path_from_experiment_index(1);
%% delete pipeline files before rerunning pipeline
% 
% todeletefiles = {'automatic_checks_complete_info.mat',...
% 'automatic_checks_complete_results.txt',...
% 'automatic_checks_incoming_info.mat',...
% 'automatic_checks_incoming_results.txt',...
% 'indicatordata.mat',...
% 'indicator_log.txt',...
% 'ledregistrationimage.png',...
% 'perframefeatures_info.mat',...
% 'registered_trx.mat',...
% 'registrationdata.mat',...
% 'registrationdata.txt',...
% 'registrationimage.png',...
% 'sexclassifier.mat',...
% 'sexclassifier_diagnostics.txt',...
% 'wingtracking_results.mat'};
% 
% todeletewildcard = {'ANALYSIS*',...
%     'ctrax_results_movie_*.mp4'};
% 
% todeletedir = {'/perframe'};
% 
% for i = 1:numel(folder_path_from_experiment_index)
%     expdir = folder_path_from_experiment_index{i};
%     % delete files
%     for j = 1:numel(todeletefiles)
%         filestr = fullfile(expdir,todeletefiles{j});
%         if exist(filestr,'file')
%             delete(filestr)
%             if exist(filestr,'file')
%                 fprintf('not deleted %s\n',filestr)
%                 return
%             end
%         end
%     end
%     %remove dirs
%     for j = 1:numel(todeletedir)
%         dirstr = fullfile(expdir,todeletedir{j});
%         
%         if exist(dirstr,'dir')
%             filelist = dir(dirstr);
%             for k = 3:numel(filelist)
%                 filestr = fullfile(expdir,todeletedir{j},filelist(k).name);
%                 
%                 if exist(filestr,'file')
%                     delete(filestr)
%                     if exist(filestr,'file')
%                         fprintf('not deleted %s\n',filestr)
%                         return
%                     end
%                 end
%             end
%             
%             rmdir(dirstr);
%             if exist(dirstr,'file')
%                 fprintf('not deleted %s\n',dirstr)
%                 return
%             end
%         end
%     end
%     % delete files with variable names
%     for j = 1:numel(todeletewildcard)
%         filelist = dir(fullfile(expdir,todeletewildcard{j}));
%         for k = 1:numel(filelist)
%             filestr = fullfile(expdir,filelist(k).name);
%             if exist(filestr,'file')
%                 delete(filestr)
%                 if exist(filestr,'file')
%                     fprintf('not deleted %s\n',filestr)
%                     return
%                 end
%             end
%         end
%     end
%     
% end


%% run analysis
% 
goldblum_analyze_experiment_folders(folder_path_from_experiment_index, settings_folder_path, lab_head_last_name, ...
                           do_force_analysis, do_use_bqueue, do_actually_submit_jobs, analysis_parameters) ;