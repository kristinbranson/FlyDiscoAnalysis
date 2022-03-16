% This isn't a proper test b/c it doesn't check whether anything worked

% Where does this script live?
this_script_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_script_path) ;
fly_disco_analysis_folder_path = fileparts(this_folder_path) ;
fly_disco_folder_path = fileparts(fly_disco_analysis_folder_path) ;
%settings_folder_path = fullfile(fly_disco_analysis_folder_path, 'settings') ;
read_only_experiments_folder_path = fullfile(fly_disco_folder_path, 'example-experiments', 'single-robie-experiment-2022-03-06-with-links-read-only') ;
working_experiments_folder_path = fullfile(fly_disco_folder_path, 'example-experiments', 'single-robie-experiment-2022-03-06-with-links') ;

% Delete the destination folder
if exist(working_experiments_folder_path, 'file') ,
    return_code = system_from_list_with_error_handling({'rm', '-rf', working_experiments_folder_path}) ;
end

% Recopy the test folder from the template
fprintf('Resetting working experiments folder...\n') ;
reset_experiment_working_copies(working_experiments_folder_path, read_only_experiments_folder_path) ;

% Find the experiments
folder_path_from_experiment_index = find_experiment_folders(working_experiments_folder_path) ;


                                
                                
                                
                                
                                
                                
                                
cluster_billing_account_name = 'branson';
do_use_bqueue = false ;    
do_actually_submit_jobs = false ;  

settings_folder_path = '/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/settings' ;
analysis_protocol = '20211014_flybubbleRed_LED';

analysis_parameters = {'analysis_protocol',analysis_protocol, ... 
    'doautomaticchecksincoming','force',...
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
    'doautomaticcheckscomplete','force'};

% folder_path_from_experiment_index = {'/groups/branson/home/robiea/Projects_data/JAABA/Data_FlyBubble/FlyTracker/cx_GMR_SS00030_CsChr_RigC_20150826T144616',...
% '/groups/branson/home/robiea/Projects_data/JAABA/Data_FlyBubble/FlyTracker/cx_GMR_SS00038_CsChr_RigB_20150729T150617',...
% '/groups/branson/home/robiea/Projects_data/JAABA/Data_FlyBubble/FlyTracker/cx_GMR_SS00168_CsChr_RigD_20150909T111218'};

% Run the script under test
fprintf('Running goldblum_analyze_experiment_folders...\n') ;
goldblum_analyze_experiment_folders(folder_path_from_experiment_index, settings_folder_path, cluster_billing_account_name, ...
                                    do_use_bqueue, do_actually_submit_jobs, analysis_parameters) ;                             
                                