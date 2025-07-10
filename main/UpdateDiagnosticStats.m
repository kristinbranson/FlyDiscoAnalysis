function UpdateDiagnosticStats(varargin)

[analysis_protocol,settingsdir,datalocparamsfilestr,...
  min_days_prev,max_days_prev] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'min_days_prev',0,'max_days_prev',28);


%% data locations

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
timestamp = datestr(now,30);
outdir = fullfile(dataloc_params.pBDPGAL4Ustatsdir,timestamp);

if ~exist(dataloc_params.pBDPGAL4Ustatsdir,'file'),
  mkdir(dataloc_params.pBDPGAL4Ustatsdir);
end

if ~exist(outdir,'file')
  mkdir(outdir);
end

%% allowed dates

d = now;
format = 'yyyymmddTHHMMSS';
maxdatenum = d - min_days_prev;
mindatenum = d - max_days_prev;
mindatestr = datestr(mindatenum,format);
maxdatestr = datestr(maxdatenum,format);

FlyBowlCombineExperiments(rootdir,outdir,...
  'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol,...
  'datalocparamsfilestr',datalocparamsfilestr,...
  'daterange',{mindatestr,maxdatestr},...
  'linename','pBDPGAL4U');