function trx = FlyBubbleComputePerFrameFeatures(expdir,varargin)

version = '0.2';

[analysis_protocol,settingsdir,datalocparamsfilestr,forcecompute,perframefns,DEBUG] = ...
  myparse(varargin,...
  'analysis_protocol','20150428_flybubble_centralcomplex',...
  'settingsdir','/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'forcecompute',true,...
  'perframefns', {}, ... % added CSC 20110321: optionally specify to-be-computed frames as parameter, reads from perframefnsfile (as before) otherwise
  'DEBUG',false...
	);

if ischar(forcecompute),
  forcecompute = str2double(forcecompute) ~= 0;
end

%% Init trx
fprintf('Initializing trx...\n');

trx = Trx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir,...
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
if forcecompute,
  WINGTRACK_PERFRAMEFILES = {'nwingsdetected' 'wing_areal' 'wing_arear' 'wing_trough_angle'};
  % AL 20131016: Blow away all preexisting (says wingtracking) to account for obsolete perframefns  
  perframefns_rm = setdiff(perframefns_preexist,WINGTRACK_PERFRAMEFILES);
  for i = 1:numel(perframefns_rm)
    pfftmp = fullfile(expdir,trx.dataloc_params.perframedir,[perframefns_rm{i} '.mat']);
    fprintf('Deleting per-frame data file %s\n',perframefns_rm{i});
    delete(pfftmp);
  end
end

%% Load trx
fprintf('Loading trajectories for %s...\n',expdir);

trx.AddExpDir(expdir,'openmovie',false); % writes trajectory fns to perframe dir

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

fprintf(logfid,'\n\n***\nRunning FlyBubbleComputePerFrameFeatures version %s analysis_protocol %s (linked to %s) at %s\n',version,analysis_protocol,real_analysis_protocol,timestamp);

%% compute per-frame features
for i = 1:nfns,
  try
    fn = perframefns{i};
    fprintf(logfid,'Computing %s...\n',fn);
    trx.(fn); 
  catch ME
    fprintf(2,'Error occurred computing fn ''%s''\n',fn);
    disp(ME.getReport());
  end    
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

