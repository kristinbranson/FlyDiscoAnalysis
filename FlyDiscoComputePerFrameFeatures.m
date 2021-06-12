function trx = FlyDiscoComputePerFrameFeatures(expdir,varargin)

version = '0.3';

p = fileparts(mfilename('fullpath'));

[analysis_protocol,settingsdir,datalocparamsfilestr,forcecompute,perframefns,DEBUG] = ...
  myparse(varargin,...
  'analysis_protocol','current_bubble',...
  'settingsdir',fullfile(p,'settings'),...
  'datalocparamsfilestr','dataloc_params.txt',...
  'forcecompute',true,...
  'perframefns', {}, ... % added CSC 20110321: optionally specify to-be-computed frames as parameter, reads from perframefnsfile (as before) otherwise
  'DEBUG',false...
	);

if ischar(forcecompute),
  forcecompute = str2double(forcecompute) ~= 0;
end

% datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
% dataloc_params = ReadParams(datalocparamsfile);

% %% 
% logger = PipelineLogger(expdir,mfilename(),dataloc_params,'perframefeature_logfilestr',...
%   settingsdir,analysis_protocol,'versionstr',version,'debug',DEBUG);

%% Init trx
fprintf('Initializing trx...\n');
trx = FBATrx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir,...
  'datalocparamsfilestr',datalocparamsfilestr,'DEBUG',DEBUG);

%% Cleanup/log existing perframefns
perframefnsfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.perframefnsfilestr);
if isempty(perframefns),
  perframefns = importdata(perframefnsfile);
end
nfns = numel(perframefns);

% what files exist already
tmp = dir(fullfile(expdir,trx.dataloc_params.perframedir,'*.mat'));
perframefns_preexist = regexprep({tmp.name},'\.mat$','');

% clean this data to force computation
WINGTRACK_PERFRAMEFILES = {'nwingsdetected' 'wing_trough_angle' 'wing_anglel' 'wing_angler'};
if trx.perframe_params.isflytracker && ~trx.perframe_params.fakectrax,
  WINGTRACK_PERFRAMEFILES = [WINGTRACK_PERFRAMEFILES,{'wing_lengthl' 'wing_lengthr'}];
else
  WINGTRACK_PERFRAMEFILES = [WINGTRACK_PERFRAMEFILES,{'wing_areal' 'wing_arear'}];
end
if forcecompute,
  
  % AL 20131016: Blow away all preexisting (sans wingtracking) to account for obsolete perframefns
  perframefns_rm = setdiff(perframefns_preexist,WINGTRACK_PERFRAMEFILES);
  for i = 1:numel(perframefns_rm)
    pfftmp = fullfile(expdir,trx.dataloc_params.perframedir,[perframefns_rm{i} '.mat']);
    fprintf('Deleting per-frame data file %s\n',perframefns_rm{i});
    delete(pfftmp);
  end
end

%% translate wing features from FlyTracker
wingperframefns = {};
if trx.perframe_params.isflytracker,
  wingperframefns = WINGTRACK_PERFRAMEFILES;
  allexist = exist(fullfile(expdir,trx.dataloc_params.wingtrxfilestr),'file');
  for i = 1:numel(WINGTRACK_PERFRAMEFILES),
    allexist = allexist && exist(fullfile(expdir,trx.dataloc_params.perframedir,[WINGTRACK_PERFRAMEFILES{i} '.mat']),'file');
    if ~allexist,
      break;
    end
  end
  if forcecompute || ~allexist,
    perframefns_rm = intersect(perframefns_preexist,WINGTRACK_PERFRAMEFILES);
    for i = 1:numel(perframefns_rm),
      pfftmp = fullfile(expdir,trx.dataloc_params.perframedir,[perframefns_rm{i} '.mat']);
      fprintf('Deleting wing tracking per-frame data file %s\n',perframefns_rm{i});
      delete(pfftmp);
    end
    
    fprintf('Running FlyTracker2WingTracking...\n');
    FlyTracker2WingTracking(expdir,'dataloc_params',trx.dataloc_params,...
      'analysis_protocol',analysis_protocol,...
      'datalocparamsfilestr',datalocparamsfilestr,...
      'settingsdir',settingsdir,...
      'perframe_params',trx.perframe_params);
  end
end

%% Load trx
fprintf('Loading trajectories for %s...\n',expdir);
trx.AddExpDir(expdir,'openmovie',false,'tryloadwingtrx',false); % writes trajectory fns to perframe dir

%% compute per-frame features
for i = 1:nfns,
  try
    fn = perframefns{i};
    fprintf('Computing %s...\n',fn);
    trx.(fn); 
  catch ME
    fprintf(2,'Error occurred computing fn ''%s''\n',fn);
    disp(ME.getReport());
  end
end

%% save info to a mat file

filename = fullfile(expdir,trx.dataloc_params.perframeinfomatfilestr);
fprintf('Saving debug info to file %s...\n',filename);
%cpffinfo = logger.runInfo;
cpffinfo = struct() ;
cpffinfo.perframefns = perframefns;
cpffinfo.forcecompute = forcecompute;
cpffinfo.perframefns_preexist = perframefns_preexist;
cpffinfo.landmark_params = trx.landmark_params;
cpffinfo.perframe_params = trx.perframe_params; 
cpffinfo.wingperframefns_computed = wingperframefns;
cpffinfo.version = version;

if exist(filename,'file'),
  try %#ok<TRYNC>
    delete(filename);
  end
end
try
  save(filename,'-struct','cpffinfo');
catch ME
  warning('FlyDiscoComputePerFrameFeatures:save',...
    'Could not save information to file %s: %s',filename,getReport(ME));
end

%% 
%logger.close();
