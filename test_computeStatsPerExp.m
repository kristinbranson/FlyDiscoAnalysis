%% Test computeStatsPerExp against reference locostatsperexp
% Run the LimbBoutAnalyzer pipeline on the same expdir, call
% computeStatsPerExp, and compare to the reference file.

modpath

%% Setup - same as Debug_FlyDiscoComputeLocomotionMetrics
settingsdir = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings';
analysis_protocol = '20251009_flybubble_LED_VNC2';
expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260112_testing_locomotioncomputeperframestats/VNC2_JRC_SS57983_RigD_20230913T120134';

%% Run FlyDiscoComputeLocomotionMetrics pipeline (replicating the setup)
fprintf('Initializing trx...\n');
trx = FBATrx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir,...
    'datalocparamsfilestr','dataloc_params.txt');
trx.AddExpDir(expdir,'dooverwrite',false,'openmovie',false);

% load aptdata
aptfile = trx.dataloc_params.apttrkfilestr;
aptdata = TrkFile.load(fullfile(expdir,aptfile));

% read stage params
stageparamsfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.locomotionmetricsparamsfilestr);
stage_params = ReadParams(stageparamsfile);
legtip_landmarknums = stage_params.legtip_landmarknums;

% load tips_velmag
load(fullfile(expdir,'tips_velmag.mat'),'tips_velmag');

% load tips_pos_body
load(fullfile(expdir,'tips_pos_body.mat'),'tips_pos_body');

% groundcontact
gc_threshold_low = stage_params.gc_threshold_low;
gc_threshold_high = stage_params.gc_threshold_high;
pairs = stage_params.pairs;
minimum_bout = stage_params.minimum_bout_groundcontact;
[groundcontact] = compute_groundcontact(tips_velmag,'pairs',pairs,...
    'gc_threshold_low',gc_threshold_low,'gc_threshold_high',gc_threshold_high,...
    'minimum_bout',minimum_bout);

% walking scores and LED indicator
[~,walking_scores] = LoadScoresFromFile(trx,'scores_Walk2',1);
indicatordata = trx.getIndicatorLED(1);
digitalindicator = indicatordata.indicatordigital;

%% Initialize LimbBoutAnalyzer and run pipeline
fprintf('Initializing LimbBoutAnalyzer...\n');
loco_analyzer = LimbBoutAnalyzer(trx, aptdata, tips_pos_body, legtip_landmarknums, ...
    groundcontact, digitalindicator, walking_scores);

fprintf('Running analyzeBoutAndStimConditions...\n');
tic;
loco_analyzer.analyzeBoutAndStimConditions();
t_bout = toc;
fprintf('  analyzeBoutAndStimConditions: %.2f seconds\n', t_bout);

fprintf('Running analyzeWalkAndStimConditions...\n');
tic;
loco_analyzer.analyzeWalkAndStimConditions();
t_walk = toc;
fprintf('  analyzeWalkAndStimConditions: %.2f seconds\n', t_walk);

%% Call computeStatsPerExp (new method)
fprintf('Running computeStatsPerExp...\n');
tic;
new_stats = loco_analyzer.computeStatsPerExp();
t_stats = toc;
fprintf('  computeStatsPerExp: %.2f seconds\n', t_stats);
fprintf('\n=== Timing Summary ===\n');
fprintf('  analyzeBoutAndStimConditions: %.2f s\n', t_bout);
fprintf('  analyzeWalkAndStimConditions: %.2f s\n', t_walk);
fprintf('  computeStatsPerExp:          %.2f s\n', t_stats);
fprintf('  Total pipeline:              %.2f s\n', t_bout + t_walk + t_stats);

%% Load reference
fprintf('Loading reference locostatsperexp...\n');
ref = load(fullfile(expdir, 'locostatsperexp.mat'));
% handle case where it may be wrapped in a variable name
if isfield(ref, 'locostatsperexp')
    ref_stats = ref.locostatsperexp;
else
    ref_stats = ref;
end

%% Compare field by field
ref_fields = fieldnames(ref_stats);
new_fields = fieldnames(new_stats);

fprintf('\n=== Comparison Results ===\n');
fprintf('Reference fields: %d\n', numel(ref_fields));
fprintf('New fields: %d\n', numel(new_fields));

% Fields in reference but not in new
missing_in_new = setdiff(ref_fields, new_fields);
if ~isempty(missing_in_new)
    fprintf('\nFields in REFERENCE but missing in NEW (%d):\n', numel(missing_in_new));
    for i = 1:numel(missing_in_new)
        fprintf('  %s\n', missing_in_new{i});
    end
end

% Fields in new but not in reference
extra_in_new = setdiff(new_fields, ref_fields);
if ~isempty(extra_in_new)
    fprintf('\nFields in NEW but not in REFERENCE (%d):\n', numel(extra_in_new));
    for i = 1:numel(extra_in_new)
        fprintf('  %s\n', extra_in_new{i});
    end
end

% Compare shared fields
shared_fields = intersect(ref_fields, new_fields);
fprintf('\nShared fields: %d\n', numel(shared_fields));

n_match = 0;
n_mismatch = 0;
tol = 1e-10;

for i = 1:numel(shared_fields)
    fn = shared_fields{i};
    ref_val = ref_stats.(fn);
    new_val = new_stats.(fn);

    match = true;
    details = '';

    % Compare .mean, .std, .Z
    for stat = {'mean', 'std', 'Z'}
        s = stat{1};
        if isfield(ref_val, s) && isfield(new_val, s)
            rv = ref_val.(s);
            nv = new_val.(s);
            if isnan(rv) && isnan(nv)
                % both NaN, OK
            elseif abs(rv - nv) > tol
                match = false;
                details = [details, sprintf(' %s: ref=%.6g new=%.6g', s, rv, nv)];
            end
        elseif isfield(ref_val, s) ~= isfield(new_val, s)
            match = false;
            details = [details, sprintf(' %s: field missing in one', s)];
        end
    end

    if match
        n_match = n_match + 1;
    else
        n_mismatch = n_mismatch + 1;
        fprintf('  MISMATCH: %s |%s\n', fn, details);
    end
end

fprintf('\n=== Summary ===\n');
fprintf('Matching fields: %d / %d\n', n_match, numel(shared_fields));
fprintf('Mismatched fields: %d / %d\n', n_mismatch, numel(shared_fields));
fprintf('Missing in new: %d\n', numel(missing_in_new));
fprintf('Extra in new: %d\n', numel(extra_in_new));

if n_mismatch == 0 && isempty(missing_in_new) && isempty(extra_in_new)
    fprintf('\nAll fields match!\n');
else
    fprintf('\nThere are differences to investigate.\n');
end
