function FlyBowlComputeTemperatureDiagnostics(expdir)

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