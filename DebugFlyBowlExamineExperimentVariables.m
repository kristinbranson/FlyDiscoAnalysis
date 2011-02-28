%% set up paths

if ispc,
  addpath E:\Code\JCtrax\misc;
  addpath E:\Code\JCtrax\filehandling;
  addpath('E:\Code\SAGE\MATLABInterface\Trunk\')
  settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
  addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
  settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
end

%% parameters

analysis_protocol = '20110222';
hfig = 1;
min_days_prev = 0;
max_days_prev = 7;

%% which experiments

% d = now;
% %format = 'yyyymmddTHHMMSS';
% format = 'yyyy-mm-ddTHH:MM:SS';
% maxdatenum = d - min_days_prev;
% mindatenum = d - max_days_prev;
% mindatestr = datestr(mindatenum,format);
% maxdatestr = datestr(maxdatenum,format);

params = {'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol,'hfig',hfig};

%% 

FlyBowlExamineExperimentVariables(params{:});