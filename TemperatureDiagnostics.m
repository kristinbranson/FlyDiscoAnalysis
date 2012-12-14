function diagnostics = TemperatureDiagnostics(expdir,varargin)


%% parse parameters
[analysis_protocol,settingsdir,datalocparamsfilestr,logfid] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'logfid',1);

version = '0.1';
timestamp = datestr(now,'yyyymmddTHHMMSS');
fprintf(logfid,'Running TemperatureDiagnostics version %s analysis_protocol %s at %s\n',version,analysis_protocol,timestamp);

%% read parameters

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

%% compute diagnostics

% read temperature file
temperaturefile = fullfile(expdir,dataloc_params.temperaturefilestr);
if ~exist(temperaturefile,'file'),
  fprintf(logfid,'Temperature stream file %s does not exist.\n',temperaturefile);
  stream = [];
else
  tempdata = importdata(temperaturefile,',');
  if isempty(tempdata),
    stream = [];
  else
    stream = tempdata(:,2);
  end
end

% diagnostics
diagnostics = struct;

% mean temperature
diagnostics.mean = nanmean(stream);

% max temperature
diagnostics.max = max(stream);

% maximum minus minimum temperature
diagnostics.maxdiff = diagnostics.max - min(stream);

% number of temperature readings
diagnostics.nreadings = numel(stream);
  
% standard deviation
diagnostics.std = nanstd(stream,1);

%% write to file
temperaturediagnosticsfile = fullfile(expdir,dataloc_params.temperaturediagnosticsfilestr);
fid = fopen(temperaturediagnosticsfile,'w');
if fid < 0,
  fprintf(logfid,'Could not open temperature diagnostics file %s for writing\n',temperaturediagnosticsfile);
else
  fns = fieldnames(diagnostics);
  for i = 1:numel(fns),
    fn = fns{i};
    fprintf(fid,'%s,%s\n',fn,num2str(diagnostics.(fn)));
  end
  fclose(fid);
end