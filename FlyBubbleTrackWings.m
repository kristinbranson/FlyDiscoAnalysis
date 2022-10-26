function [trx,perframedata,info,wtunits,trackdata] = FlyBubbleTrackWings(expdir,varargin)

% The default settings folder is "settings" within this folder
this_source_file_path = mfilename('fullpath') ;
fda_folder_path = fileparts(this_source_file_path) ;
default_settings_folder_path = fullfile(fda_folder_path, 'settings') ;

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
  'settingsdir',default_settings_folder_path,...
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
wtrxfilestr = fullfile(expdir,dataloc_params.wingtrxfilestr);
if isfield(dataloc_params,'wingperframefilestr')
    wpffilestr = fullfile(expdir,dataloc_params.wingperframefilestr);
else
    wpffilestr = fullfile(expdir,'wingtracking_perframedata.mat');
end
pfdir = fullfile(expdir,dataloc_params.perframedir);
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
    logger.log('Saving oriented trx and info to %s.\n',trxfilestr);
    save(trxfilestr,'-struct','tmp');
    % TO DO make this a relative softlink
    pwdprev = pwd;
    cd(expdir)
    cmd = sprintf('ln -s ./%s ./%s', dataloc_params.trxfilestr,dataloc_params.wingtrxfilestr);
    system(cmd);
    cd(pwdprev);
  catch ME
    warning('FlyBubbleTrackWings:save','Could not save file %s: %s\n',...
      trxfilestr,getReport(ME));
  end
end

%% save: other results
if ~DEBUG,
  if exist(wpffilestr,'file'),
    try %#ok<TRYNC>
      delete(wpffilestr);
    end
  end
  try
    logger.log('Saving wing tracking results to file %s.\n',wpffilestr);
    save(wpffilestr,'perframedata','wtunits','trackdata');
  catch ME
    warning('FlyBubbleTrackWings:save',...
      'Could not save wing tracking results to file %s: %s',wpffilestr,getReport(ME));
  end
end

%% save per-frame data directly into perframe folder

logger.log('Saving per-frame data...\n');

if ~exist(pfdir,'dir'),
  mkdir(pfdir);
end

fns = fieldnames(perframedata);
fns = setdiff(fns,{'istouching'});
for i = 1:numel(fns),
  fn = fns{i};
  s = struct('data',{perframedata.(fn)},'units',wtunits.(fn)); %#ok<NASGU>
  filename = fullfile(pfdir,[fn,'.mat']);
  save(filename,'-struct','s');
end

%% 
logger.close();
