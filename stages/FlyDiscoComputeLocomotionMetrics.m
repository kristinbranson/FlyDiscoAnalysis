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

% make list of special perframe features being computed
pfflist = {'nfeet_ground','CoM_stability'};
outputfiles = {trx.dataloc_params.locomotionmetricsswingstanceboutstatsfilestr,'tips_velmag.mat'};

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
    % need apt data 
    aptfile = trx.dataloc_params.apttrkfilestr;
    aptdata = TrkFile.load(fullfile(expdir,aptfile));
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
elseif ~exist('aptdata','var')
    aptfile = trx.dataloc_params.apttrkfilestr;
    aptdata = TrkFile.load(fullfile(expdir,aptfile));
else
    pTrk = aptdata.pTrk;
    [tips_pos_body] = compute_tips_pos_body(pTrk,apt_pts_4_center,apt_pt_4_theta,legtip_landmarknums);
    save(fullfile(expdir,'tips_pos_body.mat'),'tips_pos_body');
end


% run groundcontact detections
gc_threshold_low = stage_params.gc_threshold_low;
gc_threshold_high = stage_params.gc_threshold_high;
pairs = stage_params.pairs;
minimum_bout = stage_params.minimum_bout_groundcontact;
[groundcontact] = compute_groundcontact(tips_velmag,'pairs',pairs,'gc_threshold_low',gc_threshold_low,'gc_threshold_high',gc_threshold_high,'minimum_bout',minimum_bout);

%compute number of feet on the ground
data = cell(1,trx.nflies);
for fly = 1:trx.nflies
    data{fly} = sum(groundcontact{fly},1);
end
units = parseunits('unit');
save(fullfile(expdir,trx.dataloc_params.perframedir,'nfeet_ground.mat'),'data','units');

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

% filter swing stance for during walking
% digitscores: tarjectory format of 0,1 behavior labels
[~,walking_scores] = LoadScoresFromFile(trx,'scores_Walk2',1);

% digitalindicator = indicatorLED.indicatordigital;
indicatordata = trx.getIndicatorLED(1);
digitalindicator = indicatordata.indicatordigital;

% intialize class
loco_analyzer = LimbBoutAnalyzer(trx, aptdata, tips_pos_body, legtip_landmarknums, groundcontact, digitalindicator, walking_scores);

% compute locomotion metrics for swing, stance, and steps for walking
% during stim on and stim off periods
loco_analyzer.analyzeWalkingAndStimConditions();

% % pull bout features  bouts
% [bout_metrics_ON] = getConditionData(loco_analyzer,'walking_stimON_traj','metrics');
% [bout_metrics_OFF] = getConditionData(loco_analyzer,'walking_stimOFF_traj','metrics');

% If something goes wrong, throw an error.  Any values returned from this
% function are ignored.
try
    saveResults(loco_analyzer, fullfile(expdir,trx.dataloc_params.locomotionmetricsswingstanceboutstatsfilestr))
    % save(fullfile(expdir,trx.dataloc_params.locomotionmetricsswingstanceboutstatsfilestr),'bout_metrics_ON','bout_metrics_OFF')
catch ME
    warning('FlyDiscoComputeLocomotionMetrics:save',...
        'Could not save information to file %s: %s',filename,getReport(ME));
end
