%% set up path

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

  otherwise
    
    warning('Unknown computer %s. Paths may not be setup correctly',computername);
    addpath C:\Code\JCtrax\misc;
    addpath C:\Code\JCtrax\filehandling;
    addpath('C:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
end

%% parameters

experiment_params = struct;
% type of data to analyze
experiment_params.protocol = '*10*';
%experiment_params.rearing_protocol = '*8*';
experiment_params.screen_type = 'primary';
% no failures
experiment_params.checkflags = true;
% keep missing data, for now
experiment_params.removemissingdata = false;
% what dates should we analyze
experiment_params.daterange = {};
% name of ufmf diagnostics file
UFMFDiagnosticsFileStr = 'ufmf_diagnostics.txt';
% minimum fraction of experiment directories required to store stream
%ufmfstream_minfracexps = .9;
% mat file to save stats to 
savename = ['QuickStats_Stats_',datestr(now,30),'.mat'];

%% do it
quickstats_stats = AnalyzeQuickStats(...
  'experiment_params',struct2paramscell(experiment_params),...
  'UFMFDiagnosticsFileStr',UFMFDiagnosticsFileStr,...
  'savename',savename,...
  'rootdir',rootdatadir);