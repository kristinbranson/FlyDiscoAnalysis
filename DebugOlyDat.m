%% set up paths
[~,computername] = system('hostname');
computername = strtrim(computername);

switch computername,
  
  case 'bransonk-lw1',

    addpath E:\Code\SAGE\MATLABInterface\Trunk;
    addpath(genpath('E:\Code\cross_assay\matlab'));
    addpath(genpath('E:\Code\box\PostAnalysis'));
    addpath E:\Code\JCtrax\misc;
    addpath olydat_browser;
    rmSvnPath;
    rootdir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
    
  case 'bransonk-lw2',

    addpath C:\Code\SAGE\MATLABInterface\Trunk;
    addpath(genpath('C:\Code\cross_assay\trunk\matlab'));
    addpath(genpath('C:\Code\box\PostAnalysis\trunk\matlab'));
    addpath C:\Code\JCtrax\misc;
    addpath olydat_browser;
    rmSvnPath;
    rootdir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';

  case 'bransonk-desktop',
  
  addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
  addpath(genpath('/groups/branson/bransonlab/projects/olympiad/cross_assay/matlab'));
  addpath(genpath('/groups/branson/bransonlab/projects/olympiad/box/PostAnalysis'));
  rmSvnPath;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
  addpath olydat_browser;
  rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
  settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
  
  otherwise
    
    warning('Unknown computer %s. Paths may not be setup correctly',computername);
    addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
    addpath(genpath('/groups/branson/bransonlab/projects/olympiad/cross_assay/matlab'));
    addpath(genpath('/groups/branson/bransonlab/projects/olympiad/box/PostAnalysis'));
    rmSvnPath;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
    addpath olydat_browser;
    rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
    settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    
end

analysis_protocol = '20110407';

%% try the data selector

dsf = SAGE.Lab('olympiad').assay('bowl');
dataPull = @(x) pullBowlData(x,'dataset','data','rootdir',rootdir,'removemissingdata',false);
selector = OlyDat.DataSelector(dsf,dataPull);

%% load in data chosen with data selector

load('ControlData20110415to20110530.mat');

%% try the browser
browser = startBowlBrowser(data,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);