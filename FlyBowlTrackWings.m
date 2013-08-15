function [trx,perframedata,outtrxfile,info] = FlyBowlTrackWings(expdir,varargin)

[analysis_protocol,settingsdir,datalocparamsfilestr,paramsfile,logfid,DEBUG,leftovers] = ...
  myparse_nocheck(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'paramsfile','',...
  'logfid',[],...
  'debug',false...
	);


%% parameters

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

if isempty(paramsfile),
  paramsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.wingtrackingparamsfilestr);
end

%% start log

didopenlog = false;
if isempty(logfid),
  if isfield(dataloc_params,'trackwings_logfilestr') && ~DEBUG,
    logfile = fullfile(expdir,dataloc_params.automaticchecks_incoming_logfilestr);
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

%% get the real analysis protocol

real_analysis_protocol = GetRealAnalysisProtocol(analysis_protocol,settingsdir);


%% start

timestamp = datestr(now,'yyyymmddTHHMMSS');
fprintf(logfid,'\n\n***\nRunning TrackWings analysis_protocol %s (real analysis protocol %s) at %s\n',analysis_protocol,real_analysis_protocol,timestamp);


%% main function call

[trx,perframedata,outtrxfile,info] = ...
  TrackWings(expdir,...
  'paramsfile',paramsfile,...
  'moviefilestr',dataloc_params.moviefilestr,...
  'annfilestr',dataloc_params.annfilestr,...
  'trxfilestr',dataloc_params.trxfilestr,...
  'outtrxfilestr',dataloc_params.wingtrxfilestr,...
  'perframedir',dataloc_params.perframedir,...
  'debug',DEBUG,leftovers{:});

%% store extra info

info.analysis_protocol = analysis_protocol;
info.linked_analysis_protocol = real_analysis_protocol;
info.extra_parameters = leftovers;

if ~DEBUG,
  savefile = fullfile(expdir,dataloc_params.wingtrackinginfomatfilestr);
  if exist(savefile,'file'),
    try %#ok<TRYNC>
      delete(savefile);
    end
  end
  try
    fprintf('Saving wing tracking info to file %s.\n',savefile);
    save(savefile,'-struct','info');
  catch ME,
    warning('Could not save wing tracking info to file %s: %s',savefile,getReport(ME));
  end
end

%% print results to log file

fprintf(logfid,'Finished running FlyBowlTrackWings at %s.\n',datestr(now,'yyyymmddTHHMMSS'));

if didopenlog,
  fclose(logfid);
end