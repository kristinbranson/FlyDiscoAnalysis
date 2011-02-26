%% set up paths

if ispc,
  addpath E:\Code\JCtrax\misc;
  addpath E:\Code\JCtrax\filehandling;
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
end

%% parameters

analysis_protocol = '20110222';
hfig = 1;
min_days_prev = 0;
max_days_prev = 7;

%% which experiments

d = now;
format = 'yyyymmddTHHMMSS';
maxdatenum = d - min_days_prev;
mindatenum = d - max_days_prev;
mindatestr = datestr(mindatenum,format);
maxdatestr = datestr(maxdatenum,format);

params = {'analysis_protocol',analysis_protocol,'hfig',hfig,...
  'daterange',{mindatestr,maxdatestr}};

%% 

FlyBowlExamineExperimentVariables(params{:});