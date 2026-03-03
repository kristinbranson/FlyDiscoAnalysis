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
    rmSvnPath(false,false);
    rootdir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';

  case 'bransonk-desktop',
    
    addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
    addpath(genpath('/groups/branson/bransonlab/projects/olympiad/cross_assay/trunk/matlab'));
    addpath(genpath('/groups/branson/bransonlab/projects/olympiad/box/PostAnalysis'));
    rmSvnPath;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
    addpath olydat_browser;
    rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
    settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    
  case 'robiea-ww1'
    
    addpath C:\Users\robiea\Documents\Code_versioned\SAGE\MATLABInterface\Trunk;
    addpath(genpath('C:\Users\robiea\Documents\Code_versioned\cross_assay'));
    addpath(genpath('C:\Users\robiea\Documents\Code_versioned\PostAnalysis'));
    addpath C:\Users\robiea\Documents\Code_versioned\JCtrax\misc;
    addpath C:\Users\robiea\Documents\Code_versioned\FlyBowlAnalysis\olydat_browser;
    rmSvnPath;
    rootdir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    settingsdir = 'C:\Users\robiea\Documents\Code_versioned\FlyBowlAnalysis\settings';
    
  case 'robiea-ws'
    addpath '/groups/branson/home/robiea/Code_versioned/SAGE/MATLABInterface/Trunk';
    addpath(genpath('/groups/branson/home/robiea/Code_versioned/OlyDat/matlab'));
    addpath(genpath('/groups/branson/home/robiea/Code_versioned/OlyDat/postanalysis'));
    rmSvnPath
    addpath '/groups/branson/home/robiea/Code_versioned/JCtrax/misc';
    addpath '/groups/branson/home/robiea/Code_versioned/JCtrax/filehandling';
    addpath '/groups/branson/home/robiea/Code_versioned/FlyBowlAnalysis/olydat_browser';
    rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
    settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyBowlAnalysis/settings';
    
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

analysis_protocol = '20111005';

%% try the data selector

dsf = SAGE.Lab('olympiad').assay('bowl');
dataPull = @(x) pullBowlData(x,'dataset','score','rootdir',rootdir,'removemissingdata',false,...
  'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
selector = OlyDat.DataSelector(dsf,dataPull);

%% load in data chosen with data selector

load('ControlData20110315to20110530.mat');

%% or use SAGEGetBowlData

params = {'dataset','data','daterange',{'20111016T000000','20111023T000000'},...
  'line_name','pBDPGAL4U','checkflags',true,'removemissingdata',true,...
  'rootdir',rootdir,'analysis_protocol',analysis_protocol,'settingsdir',settingsdir};
data = SAGEGetBowlData(params{:});
data = PrepareOlyDatData(data);

%% try the browser

browser = startBowlBrowser(data,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);