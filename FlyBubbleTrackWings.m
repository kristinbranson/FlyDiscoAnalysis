function [trx,perframedata,info,wtunits,trackdata] = FlyBubbleTrackWings(expdir,varargin)

[analysis_protocol,settingsdir,datalocparamsfilestr,coparamsfile,wtparamsfile,logfid,DEBUG,leftovers] = ...
  myparse_nocheck(varargin,...
  'analysis_protocol','current_bubble',...
  'settingsdir','/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...,
  'coparamsfile','',...
  'wtparamsfile','',...
  'logfid',[],...
  'debug',false...
	);


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
  wtparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.wingtrackingparamsfilestr);
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

if isunix
  snapshotScript = which('repo_snapshot.sh');
  cmd = sprintf('%s -nocolor -brief %s',snapshotScript,settingsdir);
  [~,settingsSnapshot] = system(cmd);  
  if isdeployed
    %todo
    codeSnapshot = 'deployed';
  else
    cmd = sprintf('%s -nocolor -brief %s',snapshotScript,fileparts(mfilename('fullpath')));
    [~,codeSnapshot] = system(cmd);
  end
else
  settingsSnapshot = 'No snapshot, not running on unix.';
  codeSnapshot = settingsSnapshot;
end

%% start
timestamp = datestr(now,'yyyymmddTHHMMSS');
fprintf(logfid,'\n\n***\nRunning TrackWings analysis_protocol %s (real analysis protocol %s) at %s\n',analysis_protocol,real_analysis_protocol,timestamp);
fprintf(logfid,'\nCode snapshot: %s\n\n',codeSnapshot);
fprintf(logfid,'\nSettings snapshot: %s\n\n',settingsSnapshot);

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
info.settings_snapshot = settingsSnapshot;
info.code_snapshot = codeSnapshot;

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