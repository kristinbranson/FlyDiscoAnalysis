function FlyBubbleJAABADetect(expdir,varargin)

version = '0.2';

[analysis_protocol,settingsdir,datalocparamsfilestr,...
  forcecompute,DEBUG] = ...
  myparse(varargin,...
  'analysis_protocol','current_bubble',...
  'settingsdir','/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'forcecompute',false,...
  'debug',false);

if isunix && settingsdir(1) ~= '/',
  warning('settingsdir path must be the global path'); %#ok<*WNTAG>
end

if ischar(forcecompute),
  forcecompute = str2double(forcecompute);
end
if ischar(DEBUG),
  DEBUG = str2double(DEBUG);
end

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
paramsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.jaabadetectparamsfilestr);
params = ReadParams(paramsfile);

%% logging
logger = PipelineLogger(expdir,mfilename(),...
        dataloc_params,'jaabadetect_logfilestr',...
        settingsdir,analysis_protocol,'versionstr',version);

jaabaPFpath = fileparts(which('StartJAABA'));
versionfile = fullfile(jaabaPFpath,'version.txt');
jaaba_version = 'unknown';
if exist(versionfile,'file'),
  versionfid = fopen(versionfile,'r');
  if versionfid > 0,
    while true,
      s = fgetl(versionfid);
      if ~ischar(s),
        break;
      end
      s = strtrim(s);
      if ~isempty(s),
        jaaba_version = s;
        break;
      end
    end
    fclose(versionfid);
  end
end

logger.log('\n\n*** JAABA version %s\n',jaaba_version);

%% main loop

% loop through the classifier params files to process in order
if ~iscell(params.classifierparamsfiles),
  params.classifierparamsfiles = {params.classifierparamsfiles};
end

classifierinfo = [];
for i = 1:numel(params.classifierparamsfiles)  
  classifierparamsfile = fullfile(settingsdir,analysis_protocol,params.classifierparamsfiles{i});
  classifierinfocurr = JAABADetect(expdir,...
    'jablistfile',classifierparamsfile,...
    'forcecompute',forcecompute,...
    'debug',DEBUG);
  %'isrelativepath',true,...
  %'fnsrelative',{'featureparamfilename'});
  classifierinfo = structappend(classifierinfo,classifierinfocurr);
end

%% save info
if ~DEBUG,
  runinfo = logger.runInfo; %#ok<NASGU>
  savefile = fullfile(expdir,dataloc_params.jaabadetectinfomatfilestr);
  if exist(savefile,'file'),
    try %#ok<TRYNC>
      delete(savefile);
    end
  end
  try
    save(savefile,'classifierinfo','runinfo');
  catch ME
    warning('Could not save jaabadetect info to file %s: %s',savefile,getReport(ME));
  end  
end

%% 
logger.close();
