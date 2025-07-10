function trx = FlyDiscoComputeHOGHOFPerFrameFeatures(expdir,varargin)

% Process the args
[analysis_protocol, settingsdir, datalocparamsfilestr] = ...
  myparse(varargin,...
          'analysis_protocol','current_bubble',...
          'settingsdir', default_settings_folder_path(),...
          'datalocparamsfilestr','dataloc_params.txt',...
          'forcecompute',true,...
          'debug',false,...
          'do_run', [] ...
          );

% Read in the data locations
analysis_protocol_folder_path = fullfile(settingsdir, analysis_protocol) ;
datalocparamsfile = fullfile(analysis_protocol_folder_path, datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

% Do HOG/HOF stuff
pwdprev = pwd() ;
cleaner = onCleanup(@()(cd(pwdprev))) ;
trxfilestr = fullfile(expdir,dataloc_params.trxfilestr);
movfilestr = fullfile(expdir,dataloc_params.moviefilestr);
spacetimefeaturesdir = fileparts(which('preparePerFrameFtrs'));
cd (spacetimefeaturesdir);
fprintf('Computing HOG/HOF per-frame features...\n');
preparePerFrameFtrs(movfilestr,trxfilestr,false,false);

