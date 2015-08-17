function [trx,perframedata,info,wtunits,trackdata] = FlyBubbleTrackWings(expdir,varargin)

[analysis_protocol,...
  settingsdir,...
  datalocparamsfilestr,...
  coparamsfile,...
  wtparamsfile,...
  logfid,...
  DEBUG,...
  leftovers] = ...
  myparse_nocheck(varargin,...
  'analysis_protocol','current_bubble',...
  'settingsdir','/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...,
  'coparamsfile','',...
  'wtparamsfile','',...
  'logfid',[],...
  'debug',false);

tfcheckbuild = ischar(DEBUG) && strcmp(DEBUG,'checkbuild');
if ischar(DEBUG)
  DEBUG = str2double(DEBUG);
end

if tfcheckbuild
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
%ifofilestr = fullfile(expdir,dataloc_params.wingtrackinginfomatfilestr);
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

%% 
logger = PipelineLogger(expdir,mfilename(),dataloc_params,'trackwings_logfilestr',...
  settingsdir,analysis_protocol,'logfid',logfid,'debug',DEBUG);

%% main function call
[trx,perframedata,info,wtunits,trackdata] = ...
  ChooseOrientationsAndTrackWings(...
  trxfilestr,movfilestr,...
  'paramsfile',coparamsfile,...
  'wingtracking_params',wtparamsfile,...
  'annfile',annfilestr,...
  'debug',DEBUG,leftovers{:});
info = structmerge(info,logger.runInfo);
info.extra_parameters = leftovers;

%% save: trxfile
if ~DEBUG
  tmp = load(trxfilestr);
  tmp.trx = trx;
  tmp.wingtrackinfo = info;
  try
    fprintf('Saving oriented trx and info to %s.\n',trxfilestr);
    save(trxfilestr,'-struct','tmp');
  catch ME
    warning('FlyBubbleTrackWings:save','Could not save file %s: %s\n',...
      trxfilestr,getReport(ME));
  end
end

%% save: other results
if ~DEBUG,
  if exist(resfilestr,'file'),
    try %#ok<TRYNC>
      delete(resfilestr);
    end
  end
  try
    fprintf('Saving wing tracking results to file %s.\n',resfilestr);
    save(resfilestr,'perframedata','wtunits','trackdata');
  catch ME
    warning('FlyBubbleTrackWings:save',...
      'Could not save wing tracking results to file %s: %s',resfilestr,getReport(ME));
  end
end

%% 
logger.close();
