function FlyBowlExtraDiagnostics(expdir,varargin)

version = '0.1';
timestamp = datestr(now,'yyyymmddTHHMMSS');

[dotemperature,dobkgd,dobias,dovideo,settingsdir,...
  datalocparamsfilestr,analysis_protocol,logfid,leftovers] = ...
  myparse_nocheck(varargin,'dotemperature',true,...
  'dobkgd',true,...
  'dobias',true,...
  'dovideo',true,...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'analysis_protocol','current',...
  'logfid',[]);
if ischar(dotemperature),
  dotemperature = str2double(dotemperature)~=0;
end
if ischar(dobkgd),
  dobkgd = str2double(dotemperature)~=0;
end
if ischar(dobias),
  dobias = str2double(dotemperature)~=0;
end
if ischar(dovideo),
  dovideo = str2double(dotemperature)~=0;
end

leftovers = [leftovers,{'settingsdir',settingsdir,...
  'datalocparamsfilestr',datalocparamsfilestr,...
  'analysis_protocol',analysis_protocol}];

%% log file


datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

didopenlog = false;
if isempty(logfid),
  
  if isfield(dataloc_params,'extradiagnostics_logfilestr'),
    logfile = fullfile(expdir,dataloc_params.extradiagnostics_logfilestr);
    logfid = fopen(logfile,'a');
    if logfid < 1,
      warning('Could not open log file %s\n',logfile);
      logfid = 1;
    else
      didopenlog = true;
    end
  else
    logfid = 1;
  end
  
end

fprintf(logfid,'\n\n***\nRunning FlyBowlExtraDiagnostics version %s analysis_protocol %s at %s\n',version,analysis_protocol,timestamp);


%% compute temperature diagnostics
if dotemperature,
  fprintf(logfid,'Computing temperature diagnostics.\n');
  try
    TemperatureDiagnostics(expdir,leftovers{:},'logfid',logfid);
  catch ME,
    fprintf(logfid,'Error computing temperature diagnostics:\n%s\n',getReport(ME));
  end
else
  fprintf(logfid,'Not computing temperature diagnostics.\n');
end

%% compute background model diagnostics
if dobkgd,
  fprintf(logfid,'Computing background model diagnostics.\n');
  try
    BkgdModelDiagnostics(expdir,leftovers{:},'logfid',logfid);
  catch ME,
    fprintf(logfid,'Error computing background model diagnostics:\n%s\n',getReport(ME));
  end
else
  fprintf(logfid,'Not computing background model diagnostics.\n');
end

%% compute bias diagnostics
if dobias,
  fprintf(logfid,'Computing bias diagnostics.\n');
  try
    BowlBiasDiagnostics(expdir,leftovers{:},'logfid',logfid);
  catch ME,
    fprintf(logfid,'Error computing bias diagnostics:\n%s\n',getReport(ME));
  end
else
  fprintf(logfid,'Not computing bias diagnostics.\n');
end

%% compute average frame rate from first and last timestamps and number of frames
if dovideo,
  fprintf(logfid,'Computing video diagnostics.\n');
  try
    VideoDiagnostics(expdir,leftovers{:},'logfid',logfid);
  catch ME,
    fprintf(logfid,'Error computing video diagnostics:\n%s\n',getReport(ME));
  end
else
  fprintf(logfid,'Not computing video diagnostics.\n');
end

%% clean up

if isdeployed,
  close all;
end

%% close log file

fprintf(logfid,'Finished running FlyBowlExtraDiagnostics at %s.\n',datestr(now,'yyyymmddTHHMMSS'));

if didopenlog,
  fclose(logfid);
end