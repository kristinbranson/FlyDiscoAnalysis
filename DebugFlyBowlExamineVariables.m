%% set up paths

if ispc,
  addpath E:\Code\JCtrax\misc;
  addpath E:\Code\JCtrax\filehandling;
  addpath('E:\Code\SAGE\MATLABInterface\Trunk\')
  settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
  rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
  
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
  addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
  settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
  rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
end

%% parameters

analysis_protocol = '20110407';
vartype = 'experiment';
dataset = 'score';
% period = 2;
% loadcacheddata = true;
% datafilename = 'DataCurationSheets/ExperimentData_20110409to20110416.mat';
% username = 'bransonk';

%% which experiments

% d = now;
% format = 'yyyymmddTHHMMSS';
% maxdatenum = d - min_days_prev;
% mindatenum = d - max_days_prev;
% mindatestr = datestr(mindatenum,format);
% maxdatestr = datestr(maxdatenum,format);

params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol,...
  'rootdatadir',rootdatadir,...
  'dataset',dataset};
%   'period',period,...
%   'loadcacheddata',loadcacheddata,...
%   'datafilename',datafilename,...
%   'username',username};

%% 

[handles,data] = FlyBowlExamineVariables(vartype,params{:});
