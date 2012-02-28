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

  otherwise
    
    warning('Unknown computer %s. Paths may not be setup correctly',computername);
    addpath C:\Code\JCtrax\misc;
    addpath C:\Code\JCtrax\filehandling;
    addpath('C:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    
end

addpath olydat_browser;

%% parameters

analysis_protocol = 'current';
daterange = {'20111001T000000','20111101T000000'};

%% behavior data


format = 'yyyymmddTHHMMSS';
vartype = 'behavior';
dataset = 'data';

maxdatenum = datenum(daterange{2},format);
mindatenum = datenum(daterange{1},format);
period = maxdatenum-mindatenum;

params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol,...
  'rootdatadir',rootdatadir,...
  'dataset',dataset,...
  'checkflags',true,...
  'screen_type','primary',...
  'line_name','pBDPGAL4U',...
  'maxdatenum',maxdatenum,...
  'period',period,...
  'plotgroup','none'};
%   'maxdatenum',datenum(maxdatestr,format),...
%   'period',period,...
%   'loadcacheddata',loadcacheddata,...
%   'datafilename',datafilename,...
%   'username',username};
%   'plotgroup','set',... 

[handles,data] = FlyBowlExamineVariables(vartype,params{:});

%% score data

vartype = 'experiment';
dataset = 'score';
% period = 2;
% loadcacheddata = true;
% datafilename = 'DataCurationSheets/ExperimentData_20110409to20110416.mat';
% username = 'bransonk';

params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol,...
  'rootdatadir',rootdatadir,...
  'dataset',dataset,...
  'checkflags',true,...
  'screen_type','primary',...
  'line_name','pBDPGAL4U',...
  'maxdatenum',maxdatenum,...
  'period',period,...
  'plotgroup','none'};
%   'maxdatenum',datenum(maxdatestr,format),...
%   'period',period,...
%   'loadcacheddata',loadcacheddata,...
%   'datafilename',datafilename,...
%   'username',username};
%   'plotgroup','set',...

[handles,data] = FlyBowlExamineVariables(vartype,params{:});

%% olydat data

% params = {'dataset','data','daterange',daterange,...
%   'line_name','pBDPGAL4U','checkflags',true,'removemissingdata',true,...
%   'rootdir',rootdatadir,'analysis_protocol',analysis_protocol,'settingsdir',settingsdir};
% data = SAGEGetBowlData(params{:});
% data = PrepareOlyDatData(data);

tmp1 = load('BehaviorData_FlyBowl_20111001to20111101.mat');
tmp2 = load('ExperimentData_FlyBowl_20111001to20111101.mat');

fns1 = fieldnames(tmp1.rawdata);
fns2 = fieldnames(tmp2.rawdata);
[experiment_names,idx1,idx2] = intersect({tmp1.rawdata.experiment_name},{tmp2.rawdata.experiment_name});
data = tmp1.rawdata(idx1);
fns2not1 = setdiff(fns2,fns1);
for k = 1:numel(fns2not1),
  fn = fns2not1{k};
  for i = 1:numel(data),
    j = idx2(i);
    data(i).(fn) = tmp2.rawdata(j).(fn);
  end
end

olydata = PrepareOlyDatData(data);
olydatfilename = sprintf('OlyData_FlyBowl_%sto%s.mat',daterange{1}(1:8),daterange{2}(1:8));
save(olydatfilename,'olydata');

browser = startBowlBrowser(olydata,'settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
