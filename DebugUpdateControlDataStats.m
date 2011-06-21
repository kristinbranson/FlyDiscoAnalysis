[~,computername] = system('hostname');
computername = strtrim(computername);

switch computername,
  
  case 'bransonk-lw1',
    
    addpath E:\Code\JCtrax\misc;
    addpath E:\Code\JCtrax\filehandling;
    rootdir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
    
  case 'bransonk-lw2',
    
    addpath C:\Code\JCtrax\misc;
    addpath C:\Code\JCtrax\filehandling;
    rootdir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    
  case 'bransonk-desktop',

    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
    settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
    
  otherwise
    
    warning('Unknown computer %s. Paths may not be setup correctly',computername);
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
    settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';

end

%% parameters

analysis_protocol = '20110407';
% params = {'settingsdir',settingsdir,...
%   'analysis_protocol',analysis_protocol,...
%   'daterange',{},...
%   'experiment_protocol',{'*8*','*7*'},...
%   'rearing_protocol','*8*'};

params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol,...
  'daterange',{},...
  'experiment_protocol',{'*8*','*7*'},...
  'rearing_protocol','*8*',...
  'handling_protocol','*4*',...
  'manual_pf',{'U','P'},...
  'automated_pf',{'U','P'},...
  'flag_redo','0',...
  'flag_aborted','0',...
  'daterange',{'','20110601T000000'}};


%%

outdir = UpdateControlDataStats(rootdir,params{:});