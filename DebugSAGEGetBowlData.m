% DebugSAGEGetBowlData

if ispc,
  addpath E:\Code\JCtrax\misc;
  addpath E:\Code\JCtrax\filehandling;
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
end

%% parameters

min_days_prev = 0;
max_days_prev = 7;

%% queries

d = now;
format = 'yyyymmddTHHMMSS';
%format = 'yyyy-mm-ddTHH:MM:SS';
maxdatenum = d - min_days_prev;
mindatenum = d - max_days_prev;
mindatestr = datestr(mindatenum,format);
maxdatestr = datestr(maxdatenum,format);

% date range
queries = struct;
queries.daterange = {mindatestr,maxdatestr};

% check for failure flags
queries.checkflags = true;

% protocol 4-5
queries.experiment_protocol = {'EP0004.xml','EP0005.xml'};

% ufmf diagnostics summary or sexclassifier diagnostics
queries.data_type = {'sexclassifier_diagnostics_*','ufmf_diagnostics_summary_*'};

%% get data

params = struct2paramscell(queries);
data = SAGEGetBowlData(params{:});
