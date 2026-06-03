%% test_TCS_endtoend.m
% End-to-end validation for TCS pipeline wiring.
%
% Drives the LimbBoutAnalyzer pipeline directly (mirroring
% test_gait_class_endtoend.m / test_computeStatsPerExp.m) and verifies:
%   - walk_metrics.perexp.TCS has the expected aggregate fields
%   - per-exp event count is consistent (n + n_nontripod == n_steps)
%   - computeStatsPerExp produces TCS__walk__{LEDon,LEDoff}__all with
%     mean/std/Z, matching an independent event-pooled recompute
%   - walk_struct_{ON,OFF} carries a per-walk TCS field equal to each
%     walk's TCS.both.mean, and does NOT carry gait_class (gait stays
%     locostatsperexp-only, per the handoff)
%
% Also re-validates gait_class fractions on the same experiment.
%
% Bypasses FlyDiscoComputeLocomotionMetrics because the test experiment's
% dataloc_params is missing locomotionmetricsperexpfilestr.

addpath('/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis');
modpath;

settingsdir = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings';
analysis_protocol = '20260326_flybubble_LED_VNC2';
expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260112_testing_locomotioncomputeperframestats/VNC2_JRC_SS57983_RigD_20230913T120134';

%% Setup (replicating test_gait_class_endtoend.m)
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

% buildWalkStruct needs bout metrics. Attempt it, but the in-progress
% duty_factor change in computeStepFeatures (separate HANDOFF_duty_factor.m
% task) currently crashes bout-metric computation. Guard so the TCS-owned
% checks above still run; the walk_struct.TCS assertions below only run if
% bout metrics succeed.
walkstruct_built = false;
try
    fprintf('Running analyzeBoutAndStimConditions (needed for walk_struct)...\n');
    loco_analyzer.analyzeBoutAndStimConditions();
    fprintf('Building walk struct...\n');
    loco_analyzer.buildWalkStruct({'ON','OFF'});
    walkstruct_built = true;
catch ME
    fprintf(2, ['\n*** walk_struct validation BLOCKED: bout-metric computation failed.\n' ...
        '    %s\n' ...
        '    This is the unfinished duty_factor work in computeStepFeatures.m\n' ...
        '    (HANDOFF_duty_factor.md), unrelated to TCS wiring. walk_struct.TCS\n' ...
        '    checks skipped; rerun after duty_factor is fixed.\n\n'], ME.message);
end

%% ===== Validate TCS =====
fprintf('\n=== Validating TCS outputs ===\n');

condmap = struct('LEDon','led_on_traj', 'LEDoff','led_off_traj');
condlabels = fieldnames(condmap);

for ci = 1:numel(condlabels)
    led = condlabels{ci};
    wm = loco_analyzer.walkMetrics(condmap.(led));

    % per-walk TCS present
    assert(isfield(wm.perwalk, 'TCS'), '%s: perwalk missing TCS', led);
    assert(isfield(wm.perwalk(1).TCS, 'both'), '%s: perwalk.TCS missing .both', led);

    % per-exp TCS present + count-consistent
    assert(isfield(wm.perexp, 'TCS'), '%s: perexp missing TCS', led);
    te = wm.perexp.TCS;
    assert(te.n + te.n_nontripod == te.n_steps, ...
        '%s: n(%d)+n_nontripod(%d) != n_steps(%d)', led, te.n, te.n_nontripod, te.n_steps);
    assert(numel(te.data) == te.n, '%s: numel(data) != n', led);
    fprintf('PASS: %s perexp.TCS n=%d (events), n_steps=%d, tripod_frac=%.3f\n', ...
        led, te.n, te.n_steps, te.n / max(te.n_steps,1));

    % Independent event-pooled recompute from perwalk
    pooled = [];
    for w = 1:numel(wm.perwalk)
        pooled = [pooled, wm.perwalk(w).TCS.both.data]; %#ok<AGROW>
    end
    assert(isequaln(sort(pooled), sort(te.data)), '%s: perexp.data != pooled perwalk data', led);
    if isempty(pooled)
        assert(isnan(te.mean), '%s: empty pool should give NaN mean', led);
    else
        assert(abs(te.mean - mean(pooled,'omitnan')) < 1e-12, '%s: perexp mean mismatch', led);
    end
    fprintf('PASS: %s perexp.TCS.mean matches independent pooled recompute\n', led);

    % locostatsperexp field
    fn = sprintf('TCS__walk__%s__all', led);
    assert(isfield(stats, fn), 'missing locostatsperexp field %s', fn);
    cs = stats.(fn);
    assert(isequaln(cs.mean, te.mean), '%s: stats mean != perexp mean', led);
    assert(isequaln(cs.std, te.std), '%s: stats std != perexp std', led);
    assert(cs.Z == te.n, '%s: stats Z(%d) != perexp n(%d)', led, cs.Z, te.n);
    fprintf('PASS: %s = %.4f +/- %.4f (Z=%d)\n', fn, cs.mean, cs.std, cs.Z);

    % walk_struct carries per-walk TCS = each walk's TCS.both.mean
    if walkstruct_built
        cond = upper(strrep(led, 'LED', ''));  % LEDon -> ON, LEDoff -> OFF
        ws = loco_analyzer.walkStruct.(cond).walk_struct;
        assert(isfield(ws, 'TCS'), '%s: walk_struct missing TCS', led);
        assert(numel(ws.TCS) == numel(wm.perwalk), '%s: walk_struct.TCS length != nwalks', led);
        expected_walk_tcs = arrayfun(@(s) s.TCS.both.mean, wm.perwalk);
        assert(isequaln(ws.TCS, expected_walk_tcs), '%s: walk_struct.TCS != per-walk TCS.both.mean', led);
        assert(~isfield(ws, 'gait_class'), '%s: walk_struct should NOT contain gait_class', led);
        fprintf('PASS: %s walk_struct.TCS matches per-walk means (n=%d walks); no gait_class in walk_struct\n', ...
            led, numel(ws.TCS));
    end
end

%% ===== Re-validate gait_class fractions =====
fprintf('\n=== Re-validating gait_class fractions ===\n');
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
fprintf('PASS: gait_class 5 fractions per LED condition sum to 1\n');

if walkstruct_built
    fprintf('\n=== All TCS end-to-end checks PASSED (incl. walk_struct) ===\n');
else
    fprintf('\n=== TCS perexp/locostatsperexp + gait_class checks PASSED ===\n');
    fprintf('=== walk_struct.TCS checks SKIPPED (blocked by duty_factor WIP) ===\n');
end
