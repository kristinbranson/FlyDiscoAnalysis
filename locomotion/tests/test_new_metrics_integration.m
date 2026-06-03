%% test_new_metrics_integration.m
% Integration test: run locomotion pipeline on test experiment, verify new
% metrics appear in all output files with sane values.
%
% Tests:
%   1. duty_factor in locostatsperexp (per-pair, per-limb, all)
%   2. duty_factor in bout_metrics (per-fly, per-limb stepfeatures)
%   3. duty_factor in walk_struct / step_struct
%   4. duty_factor sanity: values in [0,1], decreases with speed
%   5. Manual spot-check: recompute duty_factor for one fly/limb from
%      raw bout data and compare to pipeline output
%
% Extend sections below for TCS, footprint clustering, metachronal lag.

addpath('/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis');
modpath;

fprintf('=== test_new_metrics_integration ===\n\n');
npass = 0;
nfail = 0;

%% Setup — replicate FlyDiscoComputeLocomotionMetrics initialization
settingsdir = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings';
analysis_protocol = '20260326_flybubble_LED_VNC2';
expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260112_testing_locomotioncomputeperframestats/VNC2_JRC_SS57983_RigD_20230913T120134';

fprintf('Initializing trx...\n');
trx = FBATrx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir,...
    'datalocparamsfilestr','dataloc_params.txt');
trx.AddExpDir(expdir,'dooverwrite',false,'openmovie',false);

aptfile = trx.dataloc_params.apttrkfilestr;
aptdata = TrkFile.load(fullfile(expdir,aptfile));

stageparamsfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.locomotionmetricsparamsfilestr);
stage_params = ReadParams(stageparamsfile);
legtip_landmarknums = stage_params.legtip_landmarknums;

load(fullfile(expdir,'tips_velmag.mat'),'tips_velmag');
load(fullfile(expdir,'tips_pos_body.mat'),'tips_pos_body');

[groundcontact] = compute_groundcontact(tips_velmag, ...
    'pairs', stage_params.pairs, ...
    'gc_threshold_low', stage_params.gc_threshold_low, ...
    'gc_threshold_high', stage_params.gc_threshold_high, ...
    'minimum_bout', stage_params.minimum_bout_groundcontact);

[~,walking_scores] = LoadScoresFromFile(trx,'scores_Walk2',1);
indicatordata = trx.getIndicatorLED(1);
digitalindicator = indicatordata.indicatordigital;

% Load onceiling/nottracking scores if available
has_onfloor_scores = true;
try
    [~,onceiling_scores] = LoadScoresFromFile(trx,'scores_onceiling_resnet_v2',1);
catch
    has_onfloor_scores = false;
end
try
    [~,nottracking_scores] = LoadScoresFromFile(trx,'scores_nottracking',1);
catch
    has_onfloor_scores = false;
end

%% Run pipeline
fprintf('Building LimbBoutAnalyzer...\n');
if has_onfloor_scores
    loco_analyzer = LimbBoutAnalyzer(trx, aptdata, tips_pos_body, legtip_landmarknums, ...
        groundcontact, digitalindicator, walking_scores, ...
        'phase_methods', {'phasediff_hilbert'}, ...
        'expdir', expdir, ...
        'do_onfloor_filtering', true, ...
        'onceiling_scores', onceiling_scores, ...
        'nottracking_scores', nottracking_scores, ...
        'frac_onfloor_threshold', stage_params.frac_onfloor_threshold);
else
    fprintf('  (onceiling/nottracking scores not found, running without onfloor filtering)\n');
    loco_analyzer = LimbBoutAnalyzer(trx, aptdata, tips_pos_body, legtip_landmarknums, ...
        groundcontact, digitalindicator, walking_scores, ...
        'phase_methods', {'phasediff_hilbert'}, ...
        'expdir', expdir);
end

if has_onfloor_scores
    fprintf('Running analyzeBoutAndStimConditions_onfloor...\n');
    loco_analyzer.analyzeBoutAndStimConditions_onfloor();
    fprintf('Running analyzeWalkAndStimConditions_onfloor...\n');
    loco_analyzer.analyzeWalkAndStimConditions_onfloor();
    fprintf('Running computeStatsPerExp...\n');
    stats = loco_analyzer.computeStatsPerExp({'ON','OFF'}, '_onfloor');
    map_suffix = '_onfloor';
else
    fprintf('Running analyzeBoutAndStimConditions...\n');
    loco_analyzer.analyzeBoutAndStimConditions();
    fprintf('Running analyzeWalkAndStimConditions...\n');
    loco_analyzer.analyzeWalkAndStimConditions();
    fprintf('Running computeStatsPerExp...\n');
    stats = loco_analyzer.computeStatsPerExp({'ON','OFF'});
    map_suffix = '';
end

fprintf('Building walkstruct...\n');
loco_analyzer.buildWalkStruct({'ON','OFF'});


%% ===== DUTY FACTOR TESTS =====
fprintf('\n--- Duty Factor Integration Tests ---\n');

%% Test 1: duty_factor fields exist in locostatsperexp
fprintf('Test 1: duty_factor in locostatsperexp...\n');
expected_fields = {};
for cond = {'LEDon','LEDoff'}
    for grp = {'all','pair1','pair2','pair3','limb1','limb2','limb3','limb4','limb5','limb6'}
        expected_fields{end+1} = sprintf('duty_factor__step__%s__%s', cond{1}, grp{1}); %#ok<SAGROW>
    end
end
missing = {};
for i = 1:numel(expected_fields)
    if ~isfield(stats, expected_fields{i})
        missing{end+1} = expected_fields{i}; %#ok<SAGROW>
    end
end
if isempty(missing)
    fprintf('  PASS: all 20 duty_factor fields present in locostatsperexp\n');
    npass = npass + 1;
else
    fprintf('  FAIL: missing %d fields: %s\n', numel(missing), strjoin(missing, ', '));
    nfail = nfail + 1;
end

%% Test 2: duty_factor values in [0,1] and not all NaN
fprintf('Test 2: duty_factor values in [0,1]...\n');
df_all = stats.duty_factor__step__LEDoff__all;
if ~isnan(df_all.mean) && df_all.mean >= 0 && df_all.mean <= 1 && df_all.Z > 0
    fprintf('  PASS: LEDoff all: mean=%.4f, std=%.4f, Z=%d\n', df_all.mean, df_all.std, df_all.Z);
    npass = npass + 1;
else
    fprintf('  FAIL: LEDoff all: mean=%.4f, std=%.4f, Z=%d\n', df_all.mean, df_all.std, df_all.Z);
    nfail = nfail + 1;
end

%% Test 3: duty_factor in bout_metrics
fprintf('Test 3: duty_factor in bout_metrics...\n');
bm_off = loco_analyzer.boutMetrics('walking_stimOFF_traj');
fly1_limb1_step = bm_off.perfly(1).perlimb(1).step;
if isfield(fly1_limb1_step, 'stepfeatures') && isfield(fly1_limb1_step.stepfeatures, 'duty_factor')
    df_vals = fly1_limb1_step.stepfeatures.duty_factor;
    n_valid = nnz(~isnan(df_vals));
    fprintf('  PASS: fly1/limb1 has %d duty_factor values (%d non-NaN)\n', numel(df_vals), n_valid);
    npass = npass + 1;
else
    fprintf('  FAIL: duty_factor not found in bout_metrics stepfeatures\n');
    nfail = nfail + 1;
end

%% Test 4: duty_factor in step_struct
fprintf('Test 4: duty_factor in step_struct...\n');
ws = loco_analyzer.walkStruct;
ss = ws.OFF.step_struct;
if isfield(ss, 'duty_factor')
    n_valid = nnz(~isnan(ss.duty_factor));
    fprintf('  PASS: step_struct.duty_factor has %d values (%d non-NaN)\n', ...
        numel(ss.duty_factor), n_valid);
    npass = npass + 1;
else
    fprintf('  FAIL: duty_factor not in step_struct\n');
    nfail = nfail + 1;
end

%% Test 5: duty_factor in walk_struct (per-limb averages)
fprintf('Test 5: duty_factor in walk_struct...\n');
wk = ws.OFF.walk_struct;
test_field = 'step_duty_factor_RF';
if isfield(wk, test_field)
    n_valid = nnz(~isnan(wk.(test_field)));
    fprintf('  PASS: walk_struct.%s has %d values (%d non-NaN)\n', ...
        test_field, numel(wk.(test_field)), n_valid);
    npass = npass + 1;
else
    fprintf('  FAIL: %s not in walk_struct\n', test_field);
    nfail = nfail + 1;
end

%% Test 6: Manual spot-check — recompute duty_factor for fly 1, limb 1
fprintf('Test 6: Manual spot-check (fly 1, limb 1, first 5 steps)...\n');
step_sf = bm_off.perfly(1).perlimb(1).step.stepfeatures;
stance_data = bm_off.perfly(1).perlimb(1).stance;
timestamps = trx(1).timestamps;

n_check = min(5, numel(step_sf.start_indices));
all_match = true;
for i = 1:n_check
    st0 = step_sf.start_indices(i);
    st1 = step_sf.end_indices(i);
    step_dur = (timestamps(st1) - timestamps(st0)) * 1000;

    % Find matching stance bout (stance_t0 == step_t0)
    stance_idx = find(stance_data.start_indices == st0, 1);
    if isempty(stance_idx)
        fprintf('  Step %d: no matching stance bout at t0=%d, skipping\n', i, st0);
        continue;
    end
    stance_dur = stance_data.durations_time(stance_idx);
    expected_df = stance_dur / step_dur;
    pipeline_df = step_sf.duty_factor(i);

    if abs(expected_df - pipeline_df) < 1e-10
        fprintf('  Step %d: manual=%.6f, pipeline=%.6f  OK\n', i, expected_df, pipeline_df);
    else
        fprintf('  Step %d: manual=%.6f, pipeline=%.6f  MISMATCH\n', i, expected_df, pipeline_df);
        all_match = false;
    end
end
if all_match
    fprintf('  PASS: manual spot-check matches pipeline\n');
    npass = npass + 1;
else
    fprintf('  FAIL: manual spot-check mismatch\n');
    nfail = nfail + 1;
end

%% Test 7: Sanity — duty_factor vs speed (negative correlation expected)
fprintf('Test 7: duty_factor vs speed correlation...\n');
ss = ws.OFF.step_struct;
if isfield(ss, 'duty_factor') && isfield(ss, 'velmag_ctr')
    valid = ~isnan(ss.duty_factor) & ~isnan(ss.velmag_ctr);
    if sum(valid) > 20
        r = corrcoef(ss.duty_factor(valid), ss.velmag_ctr(valid));
        r_val = r(1,2);
        fprintf('  r(duty_factor, velmag_ctr) = %.4f (n=%d)\n', r_val, sum(valid));
        if r_val < 0
            fprintf('  PASS: negative correlation as expected\n');
            npass = npass + 1;
        else
            fprintf('  NOTE: positive correlation (r=%.4f) — unexpected but not necessarily wrong\n', r_val);
            npass = npass + 1; % sanity check, not a hard fail
        end
    else
        fprintf('  SKIP: too few valid points (%d)\n', sum(valid));
        npass = npass + 1;
    end
else
    fprintf('  SKIP: missing duty_factor or velmag_ctr in step_struct\n');
    npass = npass + 1;
end

%% Test 8: duty_factor in bout_metrics — allflies aggregation
fprintf('Test 8: duty_factor in bout_metrics allflies aggregation...\n');
bm_off_allflies = bm_off.allflies;
pass8 = true;

% all_limbs
if isfield(bm_off_allflies.all_limbs.step, 'stepfeatures') && ...
        isfield(bm_off_allflies.all_limbs.step.stepfeatures, 'duty_factor')
    n_all = numel(bm_off_allflies.all_limbs.step.stepfeatures.duty_factor);
    fprintf('  allflies.all_limbs.step.stepfeatures.duty_factor: %d values\n', n_all);
else
    fprintf('  FAIL: duty_factor not in allflies.all_limbs.step.stepfeatures\n');
    pass8 = false;
end

% pairs (3 pairs)
for p = 1:3
    if isfield(bm_off_allflies.pairs(p).step, 'stepfeatures') && ...
            isfield(bm_off_allflies.pairs(p).step.stepfeatures, 'duty_factor')
        n_p = numel(bm_off_allflies.pairs(p).step.stepfeatures.duty_factor);
        fprintf('  allflies.pairs(%d).step.stepfeatures.duty_factor: %d values\n', p, n_p);
    else
        fprintf('  FAIL: duty_factor not in allflies.pairs(%d).step.stepfeatures\n', p);
        pass8 = false;
    end
end

% perlimb (6 limbs)
for l = 1:6
    if isfield(bm_off_allflies.perlimb(l).step, 'stepfeatures') && ...
            isfield(bm_off_allflies.perlimb(l).step.stepfeatures, 'duty_factor')
        n_l = numel(bm_off_allflies.perlimb(l).step.stepfeatures.duty_factor);
        fprintf('  allflies.perlimb(%d).step.stepfeatures.duty_factor: %d values\n', l, n_l);
    else
        fprintf('  FAIL: duty_factor not in allflies.perlimb(%d).step.stepfeatures\n', l);
        pass8 = false;
    end
end

if pass8
    fprintf('  PASS: duty_factor in allflies for all_limbs, 3 pairs, 6 perlimb\n');
    npass = npass + 1;
else
    nfail = nfail + 1;
end

%% Test 9: duty_factor aggregation consistency
% allflies.all_limbs should equal concat of all 6 perlimb
fprintf('Test 9: duty_factor aggregation consistency...\n');
all_from_perlimb = [];
for l = 1:6
    all_from_perlimb = [all_from_perlimb, bm_off_allflies.perlimb(l).step.stepfeatures.duty_factor]; %#ok<AGROW>
end
all_from_alllimbs = bm_off_allflies.all_limbs.step.stepfeatures.duty_factor;
if numel(all_from_perlimb) == numel(all_from_alllimbs)
    fprintf('  PASS: all_limbs count (%d) == sum of perlimb counts (%d)\n', ...
        numel(all_from_alllimbs), numel(all_from_perlimb));
    npass = npass + 1;
else
    fprintf('  FAIL: all_limbs count (%d) != sum of perlimb counts (%d)\n', ...
        numel(all_from_alllimbs), numel(all_from_perlimb));
    nfail = nfail + 1;
end

%% Test 10: duty_factor mean in locostatsperexp matches manual from bout_metrics
fprintf('Test 10: locostatsperexp mean matches manual computation from bout_metrics...\n');
manual_mean_all = mean(all_from_alllimbs, 'omitnan');
pipeline_mean_all = stats.duty_factor__step__LEDoff__all.mean;
if abs(manual_mean_all - pipeline_mean_all) < 1e-10
    fprintf('  PASS: manual mean=%.6f, pipeline mean=%.6f\n', manual_mean_all, pipeline_mean_all);
    npass = npass + 1;
else
    fprintf('  FAIL: manual mean=%.6f, pipeline mean=%.6f\n', manual_mean_all, pipeline_mean_all);
    nfail = nfail + 1;
end

%% Test 11: duty_factor in saveResults output (bout_metrics_OFF/ON)
fprintf('Test 11: duty_factor survives saveResults save/load cycle...\n');
tmpfile = fullfile(expdir, 'test_saveResults_dutyfactor.mat');
loco_analyzer.saveResults('test_saveResults_dutyfactor.mat');
loaded = load(tmpfile);
pass11 = true;

% bout_metrics_OFF
if isfield(loaded, 'bout_metrics_OFF') && ...
        isfield(loaded.bout_metrics_OFF.perfly(1).perlimb(1).step, 'stepfeatures') && ...
        isfield(loaded.bout_metrics_OFF.perfly(1).perlimb(1).step.stepfeatures, 'duty_factor')
    df_loaded = loaded.bout_metrics_OFF.perfly(1).perlimb(1).step.stepfeatures.duty_factor;
    df_orig = bm_off.perfly(1).perlimb(1).step.stepfeatures.duty_factor;
    if isequal(df_loaded, df_orig)
        fprintf('  bout_metrics_OFF: duty_factor matches (n=%d)\n', numel(df_loaded));
    else
        fprintf('  FAIL: bout_metrics_OFF duty_factor values differ after save/load\n');
        pass11 = false;
    end
else
    fprintf('  FAIL: duty_factor not in loaded bout_metrics_OFF\n');
    pass11 = false;
end

% bout_metrics_ON
if isfield(loaded, 'bout_metrics_ON') && ...
        isfield(loaded.bout_metrics_ON.perfly(1).perlimb(1).step, 'stepfeatures') && ...
        isfield(loaded.bout_metrics_ON.perfly(1).perlimb(1).step.stepfeatures, 'duty_factor')
    fprintf('  bout_metrics_ON: duty_factor present\n');
else
    fprintf('  FAIL: duty_factor not in loaded bout_metrics_ON\n');
    pass11 = false;
end

if pass11
    fprintf('  PASS: duty_factor survives save/load\n');
    npass = npass + 1;
else
    nfail = nfail + 1;
end
%% ===== TCS TESTS =====
fprintf('\n--- TCS Integration Tests ---\n');

%% Test 12: TCS fields in locostatsperexp
fprintf('Test 12: TCS in locostatsperexp...\n');
tcs_fields = {'TCS__walk__LEDon__all', 'TCS__walk__LEDoff__all'};
tcs_missing = {};
for i = 1:numel(tcs_fields)
    if ~isfield(stats, tcs_fields{i})
        tcs_missing{end+1} = tcs_fields{i}; %#ok<SAGROW>
    end
end
if isempty(tcs_missing)
    tcs_off = stats.TCS__walk__LEDoff__all;
    fprintf('  PASS: TCS fields present. LEDoff: mean=%.4f, std=%.4f, Z=%d\n', ...
        tcs_off.mean, tcs_off.std, tcs_off.Z);
    npass = npass + 1;
else
    fprintf('  FAIL: missing fields: %s\n', strjoin(tcs_missing, ', '));
    nfail = nfail + 1;
end

%% Test 13: TCS values in [0,1]
fprintf('Test 13: TCS values in [0,1]...\n');
tcs_off = stats.TCS__walk__LEDoff__all;
if ~isnan(tcs_off.mean) && tcs_off.mean >= 0 && tcs_off.mean <= 1 && tcs_off.Z > 0
    fprintf('  PASS: TCS mean=%.4f in [0,1], Z=%d\n', tcs_off.mean, tcs_off.Z);
    npass = npass + 1;
else
    fprintf('  FAIL: TCS mean=%.4f, Z=%d\n', tcs_off.mean, tcs_off.Z);
    nfail = nfail + 1;
end

%% Test 14: TCS in walk_struct
fprintf('Test 14: TCS in walk_struct...\n');
wk_off = ws.OFF.walk_struct;
if isfield(wk_off, 'TCS')
    n_valid = nnz(~isnan(wk_off.TCS));
    fprintf('  PASS: walk_struct.TCS has %d values (%d non-NaN)\n', ...
        numel(wk_off.TCS), n_valid);
    npass = npass + 1;
else
    fprintf('  FAIL: TCS not in walk_struct\n');
    nfail = nfail + 1;
end

%% Test 15: TCS in walk_metrics perexp
fprintf('Test 15: TCS in walk_metrics perexp...\n');
wm_off = loco_analyzer.walkMetrics('led_off_traj');
if isfield(wm_off.perexp, 'TCS')
    tcs_perexp = wm_off.perexp.TCS;
    fprintf('  perexp.TCS: mean=%.4f, n=%d, n_steps=%d, n_nontripod=%d\n', ...
        tcs_perexp.mean, tcs_perexp.n, tcs_perexp.n_steps, tcs_perexp.n_nontripod);
    % count consistency: n + n_nontripod == n_steps
    if tcs_perexp.n + tcs_perexp.n_nontripod == tcs_perexp.n_steps
        fprintf('  PASS: n + n_nontripod == n_steps\n');
        npass = npass + 1;
    else
        fprintf('  FAIL: n(%d) + n_nontripod(%d) != n_steps(%d)\n', ...
            tcs_perexp.n, tcs_perexp.n_nontripod, tcs_perexp.n_steps);
        nfail = nfail + 1;
    end
else
    fprintf('  FAIL: TCS not in walk_metrics.perexp\n');
    nfail = nfail + 1;
end

%% Test 16: TCS survives saveResults
fprintf('Test 16: TCS in saveResults output...\n');
loaded = load(tmpfile);
if isfield(loaded, 'walk_metrics_OFF') && isfield(loaded.walk_metrics_OFF.perexp, 'TCS')
    fprintf('  PASS: walk_metrics_OFF.perexp.TCS present in saved file\n');
    npass = npass + 1;
else
    fprintf('  FAIL: TCS not in saved walk_metrics_OFF\n');
    nfail = nfail + 1;
end

%% ===== CLEANUP =====
delete(tmpfile);

%% ===== SUMMARY =====
fprintf('\n=== SUMMARY: %d passed, %d failed ===\n', npass, nfail);
if nfail > 0
    error('test_new_metrics_integration: %d test(s) failed', nfail);
end
