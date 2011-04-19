%% set up path

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;

%% parameters

rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
analysis_protocol = '20110222';
datalocparamsfilestr = 'dataloc_params.txt';
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
%datetime_format = 'yyyy-mm-ddTHH:MM:SS';
datetime_format = 'yyyymmddTHHMMSS';

%% get experiment directories

data = SAGEGetBowlData('daterange',{'20110201T000000','20110316T999999'},...
  'data_type','ufmf_diagnostics_summary_stdNPxWritten',...
  'experiment_name','FlyBowl_*',...
  'dataset','score');

%% 

seconds_fliesloaded = [data.seconds_fliesloaded];
exp_datetime_str = {data.exp_datetime};
exp_datetime_num = datenum(exp_datetime_str,datetime_format);

plot(exp_datetime_num,seconds_fliesloaded,'.');
datetick('x','mmdd');
set(gca,'XTick',floor(min(exp_datetime_num))-1:ceil(max(exp_datetime_num))+1);

% bounds chosen:
min_seconds_fliesloaded = 25;
max_seconds_fliesloaded = 225;