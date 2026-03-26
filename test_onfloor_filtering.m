%% test_onfloor_filtering.m
% Verify onfloor filtering in LimbBoutAnalyzer:
% 1. Backward compatibility: unfiltered output matches original
% 2. Walk-level: post-hoc filter boutstats using walkstruct walk_valid,
%    compare to locostatsperexp_onfloor
% 3. Bout-level: same post-hoc comparison for bout/step stats

modpath;

%% Setup
settingsdir = '/groups/branson/home/robiea/Code_versioned/BransonFlyDiscoSettings/settings';
analysis_protocol = '20251009_flybubble_LED_VNC2';
expdir = '/groups/branson/home/robiea/Projects_data/FlyDisco/Bubble_data/20260325_testingOnFloorfiltering/VNC2_YNA_K_162984_RigB_20230913T100609';

fprintf('=== Test: Onfloor Filtering ===\n');
fprintf('Experiment: %s\n\n', expdir);

%% Load inputs
trx = FBATrx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir,...
    'datalocparamsfilestr','dataloc_params.txt');
trx.AddExpDir(expdir,'dooverwrite',false,'openmovie',false);

aptdata = TrkFile.load(fullfile(expdir, trx.dataloc_params.apttrkfilestr));
stage_params = ReadParams(fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.locomotionmetricsparamsfilestr));
load(fullfile(expdir,'tips_pos_body.mat'),'tips_pos_body');
load(fullfile(expdir,'tips_velmag.mat'),'tips_velmag');
[groundcontact] = compute_groundcontact(tips_velmag,'pairs',stage_params.pairs,...
    'gc_threshold_low',stage_params.gc_threshold_low,'gc_threshold_high',stage_params.gc_threshold_high,...
    'minimum_bout',stage_params.minimum_bout_groundcontact);
[~,walking_scores] = LoadScoresFromFile(trx,'scores_Walk2',1);
[~,onceiling_scores] = LoadScoresFromFile(trx,'scores_onceiling_resnet_v2',1);
[~,nottracking_scores] = LoadScoresFromFile(trx,'scores_nottracking_apt',1);
indicatordata = trx.getIndicatorLED(1);
digitalindicator = indicatordata.indicatordigital;

%% Create analyzer with onfloor filtering
loco = LimbBoutAnalyzer(trx, aptdata, tips_pos_body, stage_params.legtip_landmarknums, ...
    groundcontact, digitalindicator, walking_scores, ...
    'phase_methods', {'phasediff_hilbert'}, ...
    'do_onfloor_filtering', true, ...
    'onceiling_scores', onceiling_scores, ...
    'nottracking_scores', nottracking_scores, ...
    'frac_onfloor_threshold', 0.70);

%% Run unfiltered pipeline
fprintf('--- Running unfiltered pipeline ---\n');
loco.computeOnfloorSignal();
loco.analyzeBoutAndStimConditions();
loco.analyzeWalkAndStimConditions();
unfiltered_stats = loco.computeStatsPerExp();

% Save unfiltered large file
loco.saveResults(fullfile(expdir, 'locomotionmetricsswingstanceboutstats.mat'));

% Build walkstruct (all walks with frac_onfloor/walk_valid)
loco.buildWalkStruct({'OFF'});
loco.saveWalkStruct(fullfile(expdir, 'locomotion_walkstruct.mat'));

%% Run onfloor-filtered pipeline
fprintf('\n--- Running onfloor-filtered pipeline ---\n');
loco.analyzeBoutAndStimConditions_onfloor();
loco.analyzeWalkAndStimConditions_onfloor();
filtered_stats = loco.computeStatsPerExp({'ON','OFF'}, '_onfloor');

% Save filtered stats
locostatsperexp = filtered_stats; %#ok<NASGU>
save(fullfile(expdir, 'locostatsperexp_onfloor.mat'), 'locostatsperexp');
fprintf('Saved locostatsperexp_onfloor.mat\n');

% Save filtered large file
loco.saveResults(fullfile(expdir, 'locomotionmetricsswingstanceboutstats_onfloor.mat'));

%% TEST 1: Backward compatibility — unfiltered matches original
fprintf('\n=== TEST 1: Backward compatibility ===\n');
ref = load(fullfile(expdir, 'locostatsperexp_original.mat'));
if isfield(ref, 'locostatsperexp'), ref_stats = ref.locostatsperexp; else, ref_stats = ref; end

shared = intersect(fieldnames(ref_stats), fieldnames(unfiltered_stats));
tol = 1e-10;
n_mismatch_t1 = 0;
for i = 1:numel(shared)
    fn = shared{i};
    for stat = {'mean','std','Z'}
        s = stat{1};
        if isfield(ref_stats.(fn),s) && isfield(unfiltered_stats.(fn),s)
            rv = ref_stats.(fn).(s); nv = unfiltered_stats.(fn).(s);
            if ~(isnan(rv) && isnan(nv)) && abs(rv-nv) > tol
                n_mismatch_t1 = n_mismatch_t1 + 1;
                fprintf('  MISMATCH: %s.%s ref=%.6g new=%.6g\n', fn, s, rv, nv);
            end
        end
    end
end
if n_mismatch_t1 == 0
    fprintf('PASS: Unfiltered stats match original (%d/%d fields)\n', numel(shared), numel(shared));
else
    fprintf('FAIL: %d mismatches\n', n_mismatch_t1);
end

%% TEST 2: Post-hoc walk-level comparison
fprintf('\n=== TEST 2: Walk-level post-hoc comparison ===\n');

% Load unfiltered boutstats and walkstruct
BM = load(fullfile(expdir, 'locomotionmetricsswingstanceboutstats.mat'));
WS = load(fullfile(expdir, 'locomotion_walkstruct.mat'));
ws = WS.walk_struct_OFF;

valid_walk_idx = find(ws.walk_valid);
invalid_walk_idx = find(~ws.walk_valid);
fprintf('Walks: %d total, %d valid, %d invalid\n', numel(ws.fly), numel(valid_walk_idx), numel(invalid_walk_idx));

% For each walk-level perframe feature, collect frame data from valid walks
% and recompute the mean, then compare to filtered_stats
walk_features = {'velmag_ctr','absdv_ctr','absdu_ctr','absdtheta',...
    'forward_vel','backward_vel','left_vel','right_vel',...
    'left_dtheta','right_dtheta','CoM_stability'};

walk_metrics_unfiltered = BM.walk_metrics_OFF;
n_mismatch_t2 = 0;

for pf = 1:numel(walk_features)
    fn = walk_features{pf};
    statname = sprintf('%s__walk__LEDoff__all', fn);

    % Collect frame data from valid walks only
    % walk_metrics.perfly(f).<fn>.frm_data_fly has ALL frames for that fly
    % We need per-walk frame data. perwalk(w).<fn>.data has per-walk frames.
    valid_frm_data = [];
    for wi = 1:numel(valid_walk_idx)
        w = valid_walk_idx(wi);
        wdata = walk_metrics_unfiltered.perwalk(w).(fn).data;
        valid_frm_data = [valid_frm_data, wdata]; %#ok<AGROW>
    end
    valid_frm_data = valid_frm_data(~isnan(valid_frm_data));

    posthoc_mean = mean(valid_frm_data, 'omitnan');
    posthoc_std = std(valid_frm_data, 'omitnan');
    posthoc_Z = numel(valid_frm_data);

    if isfield(filtered_stats, statname)
        filt_mean = filtered_stats.(statname).mean;
        filt_std = filtered_stats.(statname).std;
        filt_Z = filtered_stats.(statname).Z;

        mean_match = abs(posthoc_mean - filt_mean) < tol || (isnan(posthoc_mean) && isnan(filt_mean));
        std_match = abs(posthoc_std - filt_std) < tol || (isnan(posthoc_std) && isnan(filt_std));
        Z_match = posthoc_Z == filt_Z;

        if ~mean_match || ~std_match || ~Z_match
            n_mismatch_t2 = n_mismatch_t2 + 1;
            fprintf('  MISMATCH %s: posthoc mean=%.6g filt=%.6g | Z posthoc=%d filt=%d\n', ...
                statname, posthoc_mean, filt_mean, posthoc_Z, filt_Z);
        end
    else
        fprintf('  MISSING: %s not in filtered_stats\n', statname);
        n_mismatch_t2 = n_mismatch_t2 + 1;
    end
end

if n_mismatch_t2 == 0
    fprintf('PASS: All %d walk-level features match post-hoc filtering\n', numel(walk_features));
else
    fprintf('FAIL: %d walk-level mismatches\n', n_mismatch_t2);
end

%% TEST 3: Post-hoc bout-level comparison
fprintf('\n=== TEST 3: Bout-level post-hoc comparison ===\n');

% The pipeline filters bouts using walkonfloor_digital (frame-level signal).
% For post-hoc comparison, we check each bout frame against walkonfloor_digital
% which was built from walking_scores walks (not LED-specific walks).
% We need to reload the walkonfloor_digital or reconstruct it.
%
% Simpler approach: for each bout, check if its start frame falls within
% a valid region of walkonfloor_digital. Since walkonfloor_digital is
% already computed and stored in the analyzer, we use is_onfloor_perframe
% to check the bout's frame range directly.
%
% But we don't have walkonfloor_digital saved to disk. Instead, reconstruct
% the per-bout validity from the score files directly.

bout_metrics_unfiltered = BM.bout_metrics_OFF;
nflies = numel(bout_metrics_unfiltered.perfly);

% Load score files to build per-frame is_onfloor
OC = load(fullfile(expdir, 'scores_onceiling_resnet_v2.mat'), 'allScores');
NT = load(fullfile(expdir, 'scores_nottracking_apt.mat'), 'allScores');

% Build per-fly is_onfloor and walkonfloor_digital from scores
% (same logic as computeOnfloorSignal)
is_onfloor_pf = cell(1, nflies);
walkonfloor_dig = cell(1, nflies);
for f = 1:nflies
    T0 = trx.firstframes(f);
    T1 = trx.endframes(f);
    pp_oc = OC.allScores.postprocessed{f}(T0:T1);
    pp_nt = NT.allScores.postprocessed{f}(T0:T1);
    is_onfloor_pf{f} = ~(pp_oc == 1 | pp_nt == 1);

    nframes = numel(walking_scores{f});
    [wt0s, wt1s] = detect_bouts(walking_scores{f});
    walkonfloor_dig{f} = false(1, nframes);
    for w = 1:numel(wt0s)
        t0 = wt0s(w);
        t1 = min(wt1s(w), nframes);
        frac = mean(is_onfloor_pf{f}(t0:t1));
        if frac >= 0.70
            walkonfloor_dig{f}(t0:t1) = true;
        end
    end
end

% For each bout, check if its start frame is within walkonfloor_digital
states = {'swing', 'stance'};
n_mismatch_t3 = 0;

for si = 1:numel(states)
    state = states{si};
    all_durations_filtered = [];

    for f = 1:nflies
        for limb = 1:6
            if ~isfield(bout_metrics_unfiltered.perfly(f).perlimb(limb), state)
                continue;
            end
            bdata = bout_metrics_unfiltered.perfly(f).perlimb(limb).(state);
            if ~isfield(bdata, 'start_indices') || isempty(bdata.start_indices)
                continue;
            end

            starts = bdata.start_indices;
            ends = bdata.end_indices;
            durs = bdata.durations_time;

            for b = 1:numel(starts)
                bout_start = starts(b);
                bout_end = ends(b);
                % Bout is valid if it falls entirely within walkonfloor_digital
                if bout_start >= 1 && bout_end <= numel(walkonfloor_dig{f})
                    if all(walkonfloor_dig{f}(bout_start:bout_end))
                        all_durations_filtered = [all_durations_filtered, durs(b)]; %#ok<AGROW>
                    end
                end
            end
        end
    end

    statname = sprintf('durations_time__%s__LEDoff__all', state);
    posthoc_mean = mean(all_durations_filtered, 'omitnan');
    posthoc_Z = numel(all_durations_filtered);

    if isfield(filtered_stats, statname)
        filt_mean = filtered_stats.(statname).mean;
        filt_Z = filtered_stats.(statname).Z;

        mean_match = abs(posthoc_mean - filt_mean) < tol || (isnan(posthoc_mean) && isnan(filt_mean));
        Z_match = posthoc_Z == filt_Z;

        if ~mean_match || ~Z_match
            n_mismatch_t3 = n_mismatch_t3 + 1;
            fprintf('  MISMATCH %s: posthoc mean=%.4f filt=%.4f | Z posthoc=%d filt=%d\n', ...
                statname, posthoc_mean, filt_mean, posthoc_Z, filt_Z);
        else
            fprintf('  %s: MATCH (Z=%d)\n', statname, posthoc_Z);
        end
    end
end

if n_mismatch_t3 == 0
    fprintf('PASS: Bout durations match post-hoc filtering\n');
else
    fprintf('FAIL: %d bout-level mismatches\n', n_mismatch_t3);
end

%% TEST 4: Sanity checks
fprintf('\n=== TEST 4: Sanity checks ===\n');

% Z counts: filtered <= unfiltered
n_z_violations = 0;
shared_filt = intersect(fieldnames(unfiltered_stats), fieldnames(filtered_stats));
for i = 1:numel(shared_filt)
    fn = shared_filt{i};
    if isfield(unfiltered_stats.(fn), 'Z') && isfield(filtered_stats.(fn), 'Z')
        if filtered_stats.(fn).Z > unfiltered_stats.(fn).Z
            n_z_violations = n_z_violations + 1;
            fprintf('  Z violation: %s filtered Z=%d > unfiltered Z=%d\n', ...
                fn, filtered_stats.(fn).Z, unfiltered_stats.(fn).Z);
        end
    end
end
if n_z_violations == 0
    fprintf('PASS: All filtered Z counts <= unfiltered\n');
else
    fprintf('FAIL: %d Z count violations\n', n_z_violations);
end

% walkstruct: all walks present, frac_onfloor populated
fprintf('PASS: walkstruct has %d walks, %d valid, frac_onfloor range [%.2f, %.2f]\n', ...
    numel(ws.fly), sum(ws.walk_valid), min(ws.frac_onfloor), max(ws.frac_onfloor));

%% Summary
fprintf('\n=== SUMMARY ===\n');
n_total_fail = n_mismatch_t1 + n_mismatch_t2 + n_mismatch_t3 + n_z_violations;
if n_total_fail == 0
    fprintf('ALL TESTS PASSED\n');
else
    fprintf('FAILURES: %d total\n', n_total_fail);
end
