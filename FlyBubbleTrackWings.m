function [trx,perframedata,info,wtunits,trackdata] = FlyBubbleTrackWings(expdir,varargin)

[analysis_protocol,...
  settingsdir,...
  datalocparamsfilestr,...
  coparamsfile,...
  wtparamsfile,...
  logfid,...
  DEBUG,...
  checkbuild,...
  leftovers] = ...
  myparse_nocheck(varargin,...
  'analysis_protocol','current_bubble',...
  'settingsdir','/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...,
  'coparamsfile','',...
  'wtparamsfile','',...
  'logfid',[],...
  'debug',false,...
  'checkbuild',false); 

if ischar(DEBUG)
  DEBUG = str2double(DEBUG);
end
if ischar(checkbuild)
  checkbuild = str2double(checkbuild);
end  

if checkbuild
  % just print out the build version; only works when deployed
  if ~isdeployed
    error('FlyBubbleTrackWings:deployed','Checkbuild only available when deployed.');
  end
  tmp = importdata('build.snapshot');
  cellfun(@(x)fprintf(1,'%s\n',x),tmp);
  
  trx = [];
  perframedata = [];
  info = [];
  wtunits = [];
  trackdata = [];
  return;
end

%% parameters
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

trxfilestr = fullfile(expdir,dataloc_params.trxfilestr);
movfilestr = fullfile(expdir,dataloc_params.moviefilestr);
annfilestr = fullfile(expdir,dataloc_params.annfilestr);
resfilestr = fullfile(expdir,dataloc_params.wingtrxfilestr);
ifofilestr = fullfile(expdir,dataloc_params.wingtrackinginfomatfilestr);
if isempty(coparamsfile)
  coparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.chooseorientparamsfilestr);
end
if isempty(wtparamsfile)
  wtparamsfileshort = dataloc_params.wingtrackingparamsfilestr;
  if isempty(wtparamsfileshort)
    wtparamsfile = []; % will result in use of DefaultWingTrackingParams in ChooseOrientationsAndTrackWings
  else    
    wtparamsfile = fullfile(settingsdir,analysis_protocol,wtparamsfileshort);
  end
end

%% start log
didopenlog = false;
if isempty(logfid),
  if isfield(dataloc_params,'trackwings_logfilestr') && ~DEBUG,
    logfile = fullfile(expdir,dataloc_params.automaticchecks_incoming_logfilestr);
    logfid = fopen(logfile,'a');
    if logfid < 1,
      warning('FlyBubbleTrackWings:log','Could not open log file %s\n',logfile);
      logfid = 1;
    else
      didopenlog = true;
    end
  else
    logfid = 1;
  end
end

%% get the real analysis protocol and bookkeeping info
real_analysis_protocol = GetRealAnalysisProtocol(analysis_protocol,settingsdir);
settingsSS = FlyBubbleBaR.settingssnapshot(settingsdir);
codeSS = FlyBubbleBaR.codesnapshot();

%% start
timestamp = datestr(now,'yyyymmddTHHMMSS');
fprintf(logfid,'\n\n***\nRunning TrackWings analysis_protocol %s (real analysis protocol %s) at %s\n',analysis_protocol,real_analysis_protocol,timestamp);
fprintf(logfid,'\nCode snapshot: \n');
fprintf(logfid,'%s\n',codeSS{:});
fprintf(logfid,'\nSettings snapshot: \n');
fprintf(logfid,'%s\n',settingsSS{:});

%% main function call
[trx,perframedata,info,wtunits,trackdata] = ...
  ChooseOrientationsAndTrackWings(...
  trxfilestr,movfilestr,...
  'paramsfile',coparamsfile,...
  'wingtracking_params',wtparamsfile,...
  'annfile',annfilestr,...
  'savefile',resfilestr,...
  'debug',DEBUG,leftovers{:});

%% store extra info
info.analysis_protocol = analysis_protocol;
info.linked_analysis_protocol = real_analysis_protocol;
info.extra_parameters = leftovers;
info.settings_snapshot = settingsSS;
info.code_snapshot = codeSS;

if ~DEBUG,
  if exist(ifofilestr,'file'),
    try %#ok<TRYNC>
      delete(ifofilestr);
    end
  end
  try
    fprintf('Saving wing tracking info to file %s.\n',ifofilestr);
    save(ifofilestr,'-struct','info');
  catch ME,
    warning('FlyBubbleTrackWings:info',...
      'Could not save wing tracking info to file %s: %s',ifofilestr,getReport(ME));
  end
end

%% print results to log file
fprintf(logfid,'Finished running FlyBubbleTrackWings at %s.\n',datestr(now,'yyyymmddTHHMMSS'));
if didopenlog,
  fclose(logfid);
end