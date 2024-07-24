function trx = FlyDiscoComputePerFrameFeatures(expdir,varargin)

% Declare this, even though not really used anymore
version = '0.3';

% Process the args
[analysis_protocol, settingsdir, datalocparamsfilestr, forcecompute, perframefns, DEBUG, do_run] = ...
  myparse(varargin,...
          'analysis_protocol','current_bubble',...
          'settingsdir', default_settings_folder_path(),...
          'datalocparamsfilestr','dataloc_params.txt',...
          'forcecompute',true,...
          'perframefns', {}, ... % added CSC 20110321: optionally specify to-be-computed frames as parameter, reads from perframefnsfile (as before) otherwise
          'DEBUG',false,...
          'do_run', [] ...
          );

% Init trx
fprintf('Initializing trx...\n');
trx = FBATrx('analysis_protocol',analysis_protocol, ...
             'settingsdir',settingsdir,...
             'datalocparamsfilestr',datalocparamsfilestr, ...
             'DEBUG',DEBUG);

% Cleanup/log existing perframefns
perframefnsfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.perframefnsfilestr);
if isempty(perframefns),
  perframefns = importdata(perframefnsfile);
end
nfns = numel(perframefns);
% Note that the multi-fly per-frame-feature-computation functions have been
% updated to be multiple-chamber-aware, so they only look at flies in the same
% chamber as the fly for which the PFFs are being computed.  The
% arena-relative functions did not need to be updated, because each fly's
% registered coordinates are relative to the center of the chamber it is in,
% and all the arenas are the same size.  If at some point a rig is built with
% different-sized chambers, this will have to be revisited.
% -- ALT, 2024-03-01.

% What files exist already?
tmp = dir(fullfile(expdir,trx.dataloc_params.perframedir,'*.mat'));
perframefns_preexist = regexprep({tmp.name},'\.mat$','');

% Wing-tracking per-frame features require special treatment, so enumerate
% what those are.
wing_tracking_feature_names = {'nwingsdetected' 'wing_trough_angle' 'wing_anglel' 'wing_angler'};
if trx.perframe_params.isflytracker && ~trx.perframe_params.fakectrax,
  wing_tracking_feature_names = [wing_tracking_feature_names,{'wing_lengthl' 'wing_lengthr'}];
else
  wing_tracking_feature_names = [wing_tracking_feature_names,{'wing_areal' 'wing_arear'}];
end

% If forcecompute is true, delete all preexisting PFF files 
% Apparently the individual per-frame functions are supposed to handle this,
% but apparently some don't.
if forcecompute,  
  perframefns_rm = setdiff(perframefns_preexist,wing_tracking_feature_names);
  for i = 1:numel(perframefns_rm)
    pfftmp = fullfile(expdir,trx.dataloc_params.perframedir,[perframefns_rm{i} '.mat']);
    fprintf('Deleting per-frame data file %s\n',perframefns_rm{i});
    delete(pfftmp);
  end
end

% Translate wing features from FlyTracker
wingperframefns = {};
if trx.perframe_params.isflytracker,
  wingperframefns = wing_tracking_feature_names;
  allexist = exist(fullfile(expdir,trx.dataloc_params.wingtrxfilestr),'file');
  for i = 1:numel(wing_tracking_feature_names),
    allexist = allexist && exist(fullfile(expdir,trx.dataloc_params.perframedir,[wing_tracking_feature_names{i} '.mat']),'file');
    if ~allexist,
      break;
    end
  end
  if forcecompute || ~allexist,
    perframefns_rm = intersect(perframefns_preexist,wing_tracking_feature_names);
    for i = 1:numel(perframefns_rm),
      pfftmp = fullfile(expdir,trx.dataloc_params.perframedir,[perframefns_rm{i} '.mat']);
      fprintf('Deleting wing tracking per-frame data file %s\n',perframefns_rm{i});
      delete(pfftmp);
    end
    
    % For the FlyTracker-style tracks file, if the addpflies stage is on, use the output of it.
    % If not, fall back to dataloc_params.flytrackertrackstr.
    if is_on_or_force(do_run.addpflies) ,      
      addpflies_ft_tracks_output_file_name = determine_addpflies_ft_tracks_output_file_name(trx.dataloc_params) ;
      addpflies_ft_tracks_output_file_path = fullfile(expdir, addpflies_ft_tracks_output_file_name) ;
      ft_tracks_input_file_path = addpflies_ft_tracks_output_file_path ;
    else
      ft_tracks_input_file_path = fullfile(expdir, trx.dataloc_params.flytrackertrackstr) ;
    end

    %fprintf('Running FlyTracker2WingTracking...\n');
    FlyTracker2WingTracking(expdir,...
                            'dataloc_params',trx.dataloc_params,...
                            'analysis_protocol',analysis_protocol,...
                            'settingsdir',settingsdir,...
                            'perframe_params',trx.perframe_params, ...
                            'ftfile', ft_tracks_input_file_path) ;
    %fprintf('Done running FlyTracker2WingTracking.\n');
  end
end

% Load trx
fprintf('Loading trajectories for %s...\n',expdir);
trx.AddExpDir(expdir,'openmovie',false,'tryloadwingtrx',false);  % writes trajectory fns to perframe dir

% compute per-frame features
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

% Save info to a mat file
filename = fullfile(expdir,trx.dataloc_params.perframeinfomatfilestr);
fprintf('Saving debug info to file %s...\n',filename);
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
