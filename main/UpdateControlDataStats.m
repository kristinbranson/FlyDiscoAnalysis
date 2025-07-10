function outdir = UpdateControlDataStats(rootdir,varargin)

[analysis_protocol,settingsdir,datalocparamsfilestr,...
  min_days_prev,max_days_prev,daterange,outdir,leftovers] = ...
  myparse_nocheck(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'min_days_prev',0,'max_days_prev',28,...
  'daterange',0,...
  'outdir','');

%% data locations

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
timestamp = datestr(now,30);
if isempty(outdir),
  outdir = fullfile(dataloc_params.pBDPGAL4Ustatsdir,timestamp);
end

if ~exist(dataloc_params.pBDPGAL4Ustatsdir,'file'),
  mkdir(dataloc_params.pBDPGAL4Ustatsdir);
end

didmakedir = false;
if ~exist(outdir,'file')
  mkdir(outdir);
  didmakedir = true;
end

%% allowed dates

d = now;
format = 'yyyymmddTHHMMSS';
if ~iscell(daterange) && ~isempty(daterange) && daterange == 0,
  maxdatenum = d - min_days_prev;
  mindatenum = d - max_days_prev;
  mindatestr = datestr(mindatenum,format);
  maxdatestr = datestr(maxdatenum,format);
  daterange = {mindatestr,maxdatestr};
end

success = FlyBowlCombineExperiments2(rootdir,outdir,...
  'settingsdir',settingsdir,...
  'analysis_protocol',analysis_protocol,...
  'datalocparamsfilestr',datalocparamsfilestr,...
  'daterange',daterange,...
  'line_name','pBDPGAL4U',...
  'controldatadirstr','',...
  leftovers{:});

if ~success,
  fprintf('Could not create control stats for date range %s to %s\n',daterange{:});
  if didmakedir,
    rmdir(outdir);
  end
end