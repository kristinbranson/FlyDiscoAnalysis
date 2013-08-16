function trx = FlyBowlComputePerFrameFeatures(expdir,varargin)

version = '0.2';

[analysis_protocol,settingsdir,datalocparamsfilestr,forcecompute,perframefns,DEBUG] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'forcecompute',true,...
  'perframefns', {}, ... % added CSC 20110321: optionally specify to-be-computed frames as parameter, reads from perframefnsfile (as before) otherwise
  'DEBUG',false...
	);

if ischar(forcecompute),
  forcecompute = str2double(forcecompute) ~= 0;
end

%% load the trx

fprintf('Initializing trx...\n');

trx = Trx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir,...
  'datalocparamsfilestr',datalocparamsfilestr,'DEBUG',DEBUG);

fprintf('Loading trajectories for %s...\n',expdir);

trx.AddExpDir(expdir,'openmovie',false);

%% log file

if isfield(trx.dataloc_params,'perframefeature_logfilestr') && ~DEBUG,
  logfile = fullfile(expdir,trx.dataloc_params.perframefeature_logfilestr);
  logfid = fopen(logfile,'a');
  if logfid < 1,
    warning('Could not open log file %s\n',logfile);
    logfid = 1;
  end
else
  logfid = 1;
end

timestamp = datestr(now,'yyyymmddTHHMMSS');
real_analysis_protocol = GetRealAnalysisProtocol(analysis_protocol,settingsdir);

fprintf(logfid,'\n\n***\nRunning FlyBowlComputePerFrameFeatures version %s analysis_protocol %s (linked to %s) at %s\n',version,analysis_protocol,real_analysis_protocol,timestamp);


%% compute per-frame features

perframefnsfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.perframefnsfilestr);
if isempty(perframefns),
  perframefns = importdata(perframefnsfile);
end
nfns = numel(perframefns);

% what files exist already
tmp = dir(fullfile(perframefnsfile,'*.mat'));
perframefns_preexist = regexprep({tmp.name},'\.mat$','');

% clean this data to force computation
if forcecompute,
  %deletefns = setdiff(perframefns,Trx.TrajectoryFieldNames());
  trx.CleanPerFrameData();
end

% compute each
for i = 1:nfns,
  fn = perframefns{i};
  fprintf(logfid,'Computing %s...\n',fn);
  trx.(fn); 
end

%% save info to a mat file

filename = fullfile(expdir,trx.dataloc_params.perframeinfomatfilestr);
fprintf(logfid,'Saving debug info to file %s...\n',filename);
cpffinfo = struct;
cpffinfo.perframefns = perframefns;
cpffinfo.forcecompute = forcecompute;
cpffinfo.perframefns_preexist = perframefns_preexist;
cpffinfo.version = version;
cpffinfo.analysis_protocol = analysis_protocol;
cpffinfo.linked_analysis_protocol = real_analysis_protocol;
cpffinfo.timestamp = timestamp;
cpffinfo.landmark_params = trx.landmark_params;
cpffinfo.perframe_params = trx.perframe_params; %#ok<STRNU>

if exist(filename,'file'),
  try %#ok<TRYNC>
    delete(filename);
  end
end
try
  save(filename,'-struct','cpffinfo');
catch ME,
  warning('Could not save information to file %s',filename);
end


%% close log

fprintf(logfid,'Finished running FlyBowlComputePerFrameFeatures at %s.\n',datestr(now,'yyyymmddTHHMMSS'));

if logfid > 1,
  fclose(logfid);
end

