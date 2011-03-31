function diagnostics = TemperatureDiagnostics(expdir,varargin)

%% parse parameters
[analysis_protocol,settingsdir,datalocparamsfilestr] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt');

%% read parameters

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

%% compute diagnostics

% read temperature file
temperaturefile = fullfile(expdir,dataloc_params.temperaturefilestr);
tempdata = importdata(temperaturefile,',');
stream = tempdata(:,2);

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
  warning('Could not open temperature diagnostics file %s for writing',temperaturediagnosticsfile);
else
  fns = fieldnames(diagnostics);
  for i = 1:numel(fns),
    fn = fns{i};
    fprintf(fid,'%s,%s\n',fn,num2str(diagnostics.(fn)));
  end
  fclose(fid);
end