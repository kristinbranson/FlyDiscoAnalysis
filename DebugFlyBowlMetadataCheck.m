%% set up paths

[~,computername] = system('hostname');
computername = strtrim(computername);

switch computername,
  
  case 'bransonk-lw1',
    
    addpath E:\Code\JCtrax\misc;
    addpath E:\Code\JCtrax\filehandling;
    addpath('E:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
  
  case 'bransonk-lw2',

    addpath C:\Code\JCtrax\misc;
    addpath C:\Code\JCtrax\filehandling;
    addpath('C:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    
    case 'bransonk-desktop',
        
        addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
        addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
        addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
        settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
        rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
    case 'robiea-ws',
        
        addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
        addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
        addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
        settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
        rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
        
    otherwise
    
    warning('Unknown computer %s. Paths may not be setup correctly',computername);
    addpath C:\Code\JCtrax\misc;
    addpath C:\Code\JCtrax\filehandling;
    addpath('C:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    
end

%% parameters

params = {
  'outfilename','MetaDataCheck_20110201to20111004.txt',... % set to empty to output to stdout
  'max_maxdiff_sorting_time_minutes',5,... % maximum difference in minutes between first amd last sorting time in a set
  'max_maxdiff_starvation_time_minutes',5,... % maximum difference in minutes between first amd last starvation time in a set
  'max_maxdiff_exp_datetime_minutes',5,... % maximum difference in minutes between first and last experiment start time in a set
  'min_mindt_exp_datetime_diff_sets',15,... % minimum difference in minutes between experiment start times on the same rig in different sets
  'max_maxdiff_sorting_time_perday_days',4.5,... % maximum difference in days between sorting times for first and last experiments on the same day
  'flag_aborted','0',...
  'screen_type','primary',...
  'checkflags',false,...
  'removemissingdata',false,...%};
  'daterange',{'20111001T000000'},...
  'doexpchecks',false,'dosetchecks',true,'dodatechecks',false};


%% 

[data,iserror,msgs,iserror_exp,msgs_exp,iserror_date,msgs_date] = ...
  FlyBowlMetadataCheck(params{:});