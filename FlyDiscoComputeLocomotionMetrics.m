function FlyDiscoComputeLocomotionMetrics(expdir, varargin)
% This function implements the locomotionmetrics stage.

% varargin will be a sequence of key-value pairs, including at least the keys
% 'settingsdir', 'analysis_protocol', and 'forcecompute'.  It will also
% contain any additional stage-specific key-value pairs passed as the last
% argument to FlyDiscoPipelineStage().


%version = '0.1';
% addpath /groups/branson/home/robiea/Code_versioned/locomotion_analysis/Alice
% Parse the optional arguments
[settingsdir, analysis_protocol,datalocparamsfilestr,forcecompute,~,debug] = ...
    myparse(varargin,...
    'settingsdir', default_settings_folder_path(), ...
    'analysis_protocol', 'current', ...
    'datalocparamsfilestr','dataloc_params.txt',...
    'forcecompute', false, ...
    'do_run', [],...
    'debug', []) ;

% starting stage messgage
logfid = 1;
% timestamp = datestr(now,'yyyymmddTHHMMSS');
% real_analysis_protocol = GetRealAnalysisProtocol(analysis_protocol,settingsdir);

%fprintf(logfid,'\n\n***\nRunning FlyDiscoComputeLocomotionMetrics version %s analysis_protocol %s (linked to %s) at %s\n',version,analysis_protocol,real_analysis_protocol,timestamp);

% Initialize trx class object
fprintf('Initializing trx...\n');
trx = FBATrx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir,...
    'datalocparamsfilestr',datalocparamsfilestr);
trx.AddExpDir(expdir,'dooverwrite',false,'openmovie',false);

% load aptdata
aptfile = trx.dataloc_params.apttrkfilestr;
aptdata = TrkFile.load(fullfile(expdir,aptfile));


% make list of special perframe features being computed
pfflist = {'nfeet_ground','CoM_stability','gait_class'};
outputfiles = {trx.dataloc_params.locomotionmetricsswingstanceboutstatsfilestr, ...
    trx.dataloc_params.locomotionmetricsperexpfilestr, ...
    'tips_velmag.mat','tips_pos_body.mat','groundcontact.mat'};

% if force compute is true, delete aptPFF if they exist
%should check if files exist first
% TO DO ADD  tips_velmag to list for deleting if forcecompute (don't want
% to delete now
if forcecompute,
    for i = 1:numel(pfflist) 
        curr_aptpff = fullfile(expdir,trx.dataloc_params.perframedir,[pfflist{i} '.mat']);
        if exist(curr_aptpff,'file')
            fprintf(logfid,'Deleting per-frame data file %s\n',pfflist{i});
            delete(curr_aptpff);
        end
    end
    for i = 1:numel(outputfiles)
        curr_out = fullfile(expdir,outputfiles{i});
        if exist(curr_out,'file')
            fprintf(logfid,'Deleting FlyDiscoComputeLocomotionMetrics output files: %s\n',outputfiles{i});
            delete(curr_out);
        end
    end
end


fprintf(logfid,'Computing locomotion perframe features ...\n')
% read in parameters
stageparamsfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.locomotionmetricsparamsfilestr);
stage_params = ReadParams(stageparamsfile);

% onfloor filtering is opt-in via the param file. Default false (no filtering)
% when the param is absent, for backward-compatibility with settings dirs that
% predate onceiling detection.
if isfield(stage_params,'do_onfloor_filtering')
    do_onfloor_filtering = logical(stage_params.do_onfloor_filtering);
else
    do_onfloor_filtering = false;
end
% optionally save the (large) per-fly bout/walk metrics file; default off.
if isfield(stage_params,'save_swingstancebouts')
    save_swingstancebouts = logical(stage_params.save_swingstancebouts);
else
    save_swingstancebouts = false;
end
% speed-bin edges (mm/s) for slow/med/fast conditional metrics; default [8 20].
if isfield(stage_params,'velmag_bin_edges')
    velmag_bin_edges = stage_params.velmag_bin_edges;
else
    velmag_bin_edges = [8 20];
end
% '_onfloor' suffix on map keys and output filenames when filtering is on.
key_suffix = '';
if do_onfloor_filtering
    key_suffix = '_onfloor';
end


% get leg tip velocities

% % specify leg tip points
legtip_landmarknums = stage_params.legtip_landmarknums;

% if perframe features exist load from them

% load velmag of leg tips from JAABA apt perframe features
% transform to format ground contact expects
% nfly x 6 leg tips x frames cell array
perframestr = 'apt_view1_global_velmag_';
tips_velmag =cell(1,trx.nflies);
perframeload_success = ones(1,numel(legtip_landmarknums));

for i = 1:numel(legtip_landmarknums)
    fn = [perframestr,num2str(legtip_landmarknums(i))];
    if exist(fullfile(expdir,trx.dataloc_params.perframedir,[fn,'.mat']),'file') % currently trx class can't compute apt features on the fly
        for fly = 1:trx.nflies
            tips_velmag{fly}(i,:) = trx(fly).(fn);
        end
    else
        perframeload_success(i) = 0;
    end
end

if ~all(perframeload_success) && exist(fullfile(expdir,"tips_velmag.mat"),'file') 

    % if tips_velmag already exists load them
    load(fullfile(expdir,"tips_velmag.mat"),'tips_velmag');

else
    % compute and save if doesn't already exist
    % % need apt data 
    % aptfile = trx.dataloc_params.apttrkfilestr;
    % aptdata = TrkFile.load(fullfile(expdir,aptfile));
    ts = trx.movie_timestamps{:};
    dt = diff(ts)';
    pxpermm = trx.pxpermm;
  
        % leg tip velocity
    [tips_velmag] = compute_legtipvelmag(aptdata,dt,pxpermm,legtip_landmarknums);

    save(fullfile(expdir,'tips_velmag.mat'),'tips_velmag');
end

% tips_pos_body
apt_pts_4_center = stage_params.apt_pts_4_center;
apt_pt_4_theta = stage_params.apt_pt_4_theta;
if exist(fullfile(expdir,"tips_pos_body.mat"),'file')
    % if tips positions in body reference already exists load them
    load(fullfile(expdir,"tips_pos_body.mat"),'tips_pos_body');
else
    % if ~exist('aptdata','var')
    %     aptfile = trx.dataloc_params.apttrkfilestr;
    %     aptdata = TrkFile.load(fullfile(expdir,aptfile));
    % end
    pTrk = aptdata.pTrk;
    [tips_pos_body] = compute_tips_pos_body(pTrk,apt_pts_4_center,apt_pt_4_theta,legtip_landmarknums);
    save(fullfile(expdir,'tips_pos_body.mat'),'tips_pos_body');
end


% run groundcontact detections
gc_threshold_low = stage_params.gc_threshold_low;
gc_threshold_high = stage_params.gc_threshold_high;
pairs = stage_params.pairs;
minimum_bout = stage_params.minimum_bout_groundcontact;
groundcontact_file = fullfile(expdir,'groundcontact.mat');
if exist(groundcontact_file,'file') && ~forcecompute
    load(groundcontact_file,'groundcontact');
else
    [groundcontact] = compute_groundcontact(tips_velmag,'pairs',pairs,'gc_threshold_low',gc_threshold_low,'gc_threshold_high',gc_threshold_high,'minimum_bout',minimum_bout);
    save(groundcontact_file,'groundcontact');
end

%compute number of feet on the ground
data = cell(1,trx.nflies);
for fly = 1:trx.nflies
    data{fly} = sum(groundcontact{fly},1);
end
units = parseunits('unit');
save(fullfile(expdir,trx.dataloc_params.perframedir,'nfeet_ground.mat'),'data','units');

% compute gait class per frame (1=tripod, 2=tetrapod, 3=grounded,
% 4=airborne, 5=other). See gait_pattern_constants.m.
gaitfilestr = fullfile(expdir,trx.dataloc_params.perframedir,'gait_class.mat');
if exist(gaitfilestr,'file') && ~forcecompute
    % already on disk; nothing to do (LimbBoutAnalyzer reads via trx)
else
    data = compute_gait_class(groundcontact);
    units = parseunits('unit');
    save(gaitfilestr,'data','units');
end

% compute COM took ~24 seconds to compute, sped up with parfor ~3 seconds.
% but WAY slower if have to start parallel pool. 
CoMfilestr = fullfile(expdir,trx.dataloc_params.perframedir,'CoM_stability.mat');
if exist(CoMfilestr,'file')
    % if tips positions in body reference already exists load them
    load(CoMfilestr,'data');
    CoM_stability = data;
else
    CoM_stabilty = compute_CoMstability(aptdata,legtip_landmarknums,groundcontact);
    data = CoM_stabilty;
    units = parseunits('px');
    save(fullfile(expdir,trx.dataloc_params.perframedir,'CoM_stability.mat'),'data','units');
end

fprintf(logfid,'Computing swing and stance bout metrics ...\n')

%compute swing stance bouts
% [perfly_limbboutdata] = limbSwingStanceStep(groundcontact);

% Load walking scores from stage_params
[~,walking_scores] = LoadScoresFromFile(trx, stage_params.walking_score_file, 1);

% Load onfloor filtering scores only when filtering is enabled. When enabled,
% the score files are required: LoadScoresFromFile errors if they are missing
% (deterministic - do not silently skip filtering).
if do_onfloor_filtering
    fprintf(logfid,'onfloor filtering ON: requiring onceiling + nottracking scores\n');
    [~,onceiling_scores] = LoadScoresFromFile(trx, stage_params.onceiling_score_file, 1);
    [~,nottracking_scores] = LoadScoresFromFile(trx, stage_params.nottracking_score_file, 1);
else
    fprintf(logfid,'onfloor filtering OFF: walks are not filtered for onceiling/nottracking\n');
    onceiling_scores = {};
    nottracking_scores = {};
end

% digitalindicator
indicatordata = trx.getIndicatorLED(1);
digitalindicator = indicatordata.indicatordigital;

% Initialize analyzer (onfloor filtering optional)
loco_analyzer = LimbBoutAnalyzer(trx, aptdata, tips_pos_body, legtip_landmarknums, groundcontact, digitalindicator, walking_scores, ...
    'phase_methods', {'phasediff_hilbert'}, ...
    'do_onfloor_filtering', do_onfloor_filtering, ...
    'onceiling_scores', onceiling_scores, ...
    'nottracking_scores', nottracking_scores, ...
    'frac_onfloor_threshold', stage_params.frac_onfloor_threshold);
loco_analyzer.velmag_bin_edges = velmag_bin_edges;

% compute bout + walk metrics for walking during stim on/off, using the
% onfloor-filtered or plain methods depending on do_onfloor_filtering.
if do_onfloor_filtering
    loco_analyzer.analyzeBoutAndStimConditions_onfloor();
    loco_analyzer.analyzeWalkAndStimConditions_onfloor();
else
    loco_analyzer.analyzeBoutAndStimConditions();
    loco_analyzer.analyzeWalkAndStimConditions();
end

% compute and save locostatsperexp; output filename reflects filtering mode.
loco_analyzer.computeStatsPerExp({'ON','OFF'}, key_suffix);
try
    perexpfilename = apply_mode_suffix(trx.dataloc_params.locomotionmetricsperexpfilestr, key_suffix);
    loco_analyzer.saveStatsPerExp(fullfile(expdir, perexpfilename));
catch ME
    warning('FlyDiscoComputeLocomotionMetrics:saveStatsPerExp',...
        'Could not save per-experiment stats to file %s: %s',perexpfilename,getReport(ME));
end

% optionally save the full per-fly bout/walk metrics (large file).
if save_swingstancebouts
    try
        ssbfilename = apply_mode_suffix(trx.dataloc_params.locomotionmetricsswingstanceboutstatsfilestr, key_suffix);
        loco_analyzer.saveResults(fullfile(expdir, ssbfilename), key_suffix);
    catch ME
        warning('FlyDiscoComputeLocomotionMetrics:saveResults',...
            'Could not save swingstancebouts file: %s',getReport(ME));
    end
end

end  % function FlyDiscoComputeLocomotionMetrics


function fname = apply_mode_suffix(fname, key_suffix)
% Insert key_suffix before the extension, after stripping any pre-existing
% '_onfloor' baked into the dataloc filename, so the name reflects the mode
% actually run (e.g. 'locostatsperexp_onfloor.mat' -> base 'locostatsperexp').
[d,nm,ext] = fileparts(fname);
nm = regexprep(nm,'_onfloor$','');
fname = fullfile(d,[nm key_suffix ext]);
end
