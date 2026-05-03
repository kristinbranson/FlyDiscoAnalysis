%% test_gait_class_endtoend.m
% End-to-end validation for Phase C gait_class wiring.
%
% Drives the LimbBoutAnalyzer pipeline directly (mirroring
% test_computeStatsPerExp.m setup) and verifies:
%   - perframe gait_class file path used by the fallback works
%   - walk_metrics.perexp.gait_class has the expected count fields
%   - computeStatsPerExp produces 5 gait_class fraction fields per LED
%     condition and they sum to 1
%   - nfeet_ground also flows into the per-experiment stats
%
% This bypasses FlyDiscoComputeLocomotionMetrics because the test
% experiment's dataloc_params is missing locomotionmetricsperexpfilestr.

addpath('/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis');
modpath;

settingsdir = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings';
analysis_protocol = '20251009_flybubble_LED_VNC2';
expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260112_testing_locomotioncomputeperframestats/VNC2_JRC_SS57983_RigD_20230913T120134';

% Wipe gait_class.mat from perframe dir to exercise the in-memory
% fallback path inside compute_WalkFeatures.
gait_pf = fullfile(expdir, 'perframe', 'gait_class.mat');
if exist(gait_pf, 'file')
    delete(gait_pf);
    fprintf('Deleted existing %s (testing fallback path)\n', gait_pf);
end

%% Setup (replicating test_computeStatsPerExp.m)
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

%% Build analyzer + run walk pipeline
fprintf('Building LimbBoutAnalyzer...\n');
loco_analyzer = LimbBoutAnalyzer(trx, aptdata, tips_pos_body, legtip_landmarknums, ...
    groundcontact, digitalindicator, walking_scores, ...
    'phase_methods', {'phasediff_hilbert'}, ...
    'expdir', expdir);

fprintf('Running analyzeWalkAndStimConditions...\n');
loco_analyzer.analyzeWalkAndStimConditions();

fprintf('Running computeStatsPerExp...\n');
stats = loco_analyzer.computeStatsPerExp({'ON','OFF'});

%% Validation
fprintf('\n=== Validating gait_class outputs ===\n');

% Per-walk struct should have gait_class with .data and 5 counts
wm_off = loco_analyzer.walkMetrics('led_off_traj');
assert(isfield(wm_off.perwalk, 'gait_class'), 'perwalk missing gait_class');
gc1 = wm_off.perwalk(1).gait_class;
required = {'data','n_frames','tripod_count','tetrapod_count', ...
            'grounded_count','airborne_count','other_count'};
for i = 1:numel(required)
    assert(isfield(gc1, required{i}), 'perwalk(1).gait_class missing %s', required{i});
end
fprintf('PASS: perwalk(1).gait_class has fields %s\n', strjoin(fieldnames(gc1), ', '));

% Class codes sane: integers in [1,5], counts add to n_frames
assert(all(gc1.data >= 1 & gc1.data <= 5), 'class codes out of [1,5]');
total = gc1.tripod_count + gc1.tetrapod_count + gc1.grounded_count + ...
        gc1.airborne_count + gc1.other_count;
assert(total == gc1.n_frames, 'count sum (%d) != n_frames (%d)', total, gc1.n_frames);
fprintf('PASS: walk-level counts consistent (n_frames=%d)\n', gc1.n_frames);

% Per-fly struct
assert(isfield(wm_off.perfly, 'gait_class'), 'perfly missing gait_class');
gcf = wm_off.perfly(1).gait_class;
assert(numel(gcf.frm_data_fly) == gcf.frm_n_fly, 'frm_data_fly/frm_n_fly mismatch');
fprintf('PASS: perfly(1).gait_class.frm_n_fly = %d\n', gcf.frm_n_fly);

% Per-experiment counts
assert(isfield(wm_off.perexp, 'gait_class'), 'perexp missing gait_class');
gce = wm_off.perexp.gait_class;
total_exp = gce.tripod_count_exp + gce.tetrapod_count_exp + gce.grounded_count_exp + ...
            gce.airborne_count_exp + gce.other_count_exp;
assert(total_exp == gce.frm_n_exp, 'per-exp count sum (%d) != frm_n_exp (%d)', ...
       total_exp, gce.frm_n_exp);
fprintf('PASS: per-exp counts consistent (frm_n_exp=%d)\n', gce.frm_n_exp);

% Per-experiment fraction fields
labels = {'tripod','tetrapod','grounded','airborne','other'};
for cond = {'LEDoff','LEDon'}
    tot_frac = 0;
    fprintf('  %s: ', cond{1});
    for i = 1:numel(labels)
        fn = sprintf('gait_class__walk__%s__%s_frac', cond{1}, labels{i});
        assert(isfield(stats, fn), 'missing field %s', fn);
        tot_frac = tot_frac + stats.(fn).frac;
        fprintf('%s=%.3f  ', labels{i}, stats.(fn).frac);
    end
    fprintf('| sum=%.6f\n', tot_frac);
    assert(abs(tot_frac - 1) < 1e-9, '%s fractions sum to %g, expected 1', cond{1}, tot_frac);
end
fprintf('PASS: 5 fraction fields per LED condition; sums = 1\n');

% nfeet_ground per-experiment field present
nfeet_fn = 'nfeet_ground__walk__LEDoff__all';
assert(isfield(stats, nfeet_fn), 'expected field %s missing', nfeet_fn);
fprintf('PASS: %s.mean = %.4f, .std = %.4f, .Z = %d\n', ...
        nfeet_fn, stats.(nfeet_fn).mean, stats.(nfeet_fn).std, stats.(nfeet_fn).Z);

fprintf('\n=== All Phase C end-to-end checks PASSED ===\n');
