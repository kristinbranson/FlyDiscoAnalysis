%% test_TCS_aggregation.m
% Synthetic-data unit test for the TCS aggregation layer:
%   computePerFlyTCS and computePerExpTCS.
%
% Exercises only concatenation + linear mean/std/n + count summation. The
% underlying per-walk TCS computation is covered by test_compute_TCS.m, and
% combineTCSMetrics + field naming are covered end-to-end by
% test_TCS_endtoend.m.
%
% Coverage:
%   1. computePerFlyTCS: concat data/event_speeds across walks, sum
%      n_steps/n_nontripod, linear mean/std/n.
%   2. computePerFlyTCS: walks where some have NO tripod events (empty data
%      but n_steps>0) -> events from the non-empty walks only.
%   3. computePerFlyTCS: fly with NO walks (0-element struct) -> NaN mean,
%      n=0, n_steps=0.
%   4. computePerExpTCS: concat across walks (= across flies).
%   5. computePerExpTCS: experiment with no tripod events anywhere ->
%      NaN mean, n=0.

addpath('/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis');
modpath;

fprintf('=== test_TCS_aggregation ===\n\n');

%% --- Test 1: computePerFlyTCS over two walks with events -------------------
% walk1: data=[0.8 0.9], walk2: data=[0.7]
wfs = make_walkstruct({[0.8 0.9], [0.7]}, ...
                      {[10 12], [14]}, ...   % event_speeds
                      [3 2], [1 1]);         % n_steps, n_nontripod per walk

perfly = computePerFlyTCS(wfs);

expected_data = [0.8 0.9 0.7];
assert(isequal(perfly.data, expected_data), 'Test 1: data concat mismatch');
assert(isequal(perfly.event_speeds, [10 12 14]), 'Test 1: event_speeds concat mismatch');
assert(abs(perfly.mean - mean(expected_data)) < 1e-12, 'Test 1: mean mismatch');
assert(abs(perfly.std - std(expected_data)) < 1e-12, 'Test 1: std mismatch');
assert(perfly.n == 3, 'Test 1: n expected 3, got %d', perfly.n);
assert(perfly.n_steps == 5, 'Test 1: n_steps expected 5, got %d', perfly.n_steps);
assert(perfly.n_nontripod == 2, 'Test 1: n_nontripod expected 2, got %d', perfly.n_nontripod);
fprintf('Test 1 PASSED: computePerFlyTCS concat + linear stats + count sums\n');

%% --- Test 2: computePerFlyTCS with a walk that has no tripod events --------
% walk1: 2 events; walk2: 0 events (all candidate steps non-tripod);
% walk3: 1 event. Empty walk must still contribute its n_steps/n_nontripod.
wfs2 = make_walkstruct({[0.8 0.9], [], [0.6]}, ...
                       {[10 12], [], [8]}, ...
                       [3 4 2], [1 4 1]);

perfly2 = computePerFlyTCS(wfs2);
expected2 = [0.8 0.9 0.6];
assert(isequal(perfly2.data, expected2), 'Test 2: data should skip empty walk');
assert(isequal(perfly2.event_speeds, [10 12 8]), 'Test 2: event_speeds mismatch');
assert(abs(perfly2.mean - mean(expected2)) < 1e-12, 'Test 2: mean mismatch');
assert(perfly2.n == 3, 'Test 2: n expected 3, got %d', perfly2.n);
assert(perfly2.n_steps == 9, 'Test 2: n_steps expected 9, got %d', perfly2.n_steps);
assert(perfly2.n_nontripod == 6, 'Test 2: n_nontripod expected 6, got %d', perfly2.n_nontripod);
fprintf('Test 2 PASSED: walk with no tripod events contributes counts only\n');

%% --- Test 3: computePerFlyTCS for a fly with NO walks ----------------------
% Degenerate 0-element walkfeaturestruct (the pipeline skips such flies, but
% the aggregator must not error).
wfs_nowalks = make_walkstruct({}, {}, [], []);
assert(numel(wfs_nowalks) == 0, 'Test 3 setup: expected 0-element struct');
perfly3 = computePerFlyTCS(wfs_nowalks);
assert(isempty(perfly3.data), 'Test 3: data should be empty');
assert(isnan(perfly3.mean), 'Test 3: mean should be NaN');
assert(perfly3.n == 0, 'Test 3: n should be 0');
assert(perfly3.n_steps == 0, 'Test 3: n_steps should be 0');
assert(perfly3.n_nontripod == 0, 'Test 3: n_nontripod should be 0');
fprintf('Test 3 PASSED: fly with no walks -> empty/NaN, counts 0\n');

%% --- Test 4: computePerExpTCS over walks from two flies --------------------
% fly1 walks: [0.8 0.9], [0.7];  fly2 walks: [0.6 0.5]
% perwalk_metrics is the flat walk array across all flies.
perwalk = make_walkstruct({[0.8 0.9], [0.7], [0.6 0.5]}, ...
                          {[10 12], [14], [8 9]}, ...
                          [3 2 2], [1 1 0]);
perfly_dummy = struct('TCS', {[]});  % unused by computePerExpTCS

perexp = computePerExpTCS(perwalk, perfly_dummy);

expected_all = [0.8 0.9 0.7 0.6 0.5];
assert(isequal(perexp.data, expected_all), 'Test 4: data concat mismatch');
assert(isequal(perexp.event_speeds, [10 12 14 8 9]), 'Test 4: event_speeds concat mismatch');
assert(abs(perexp.mean - mean(expected_all)) < 1e-12, 'Test 4: mean mismatch');
assert(abs(perexp.std - std(expected_all)) < 1e-12, 'Test 4: std mismatch');
assert(perexp.n == 5, 'Test 4: n expected 5, got %d', perexp.n);
assert(perexp.n_steps == 7, 'Test 4: n_steps expected 7, got %d', perexp.n_steps);
assert(perexp.n_nontripod == 2, 'Test 4: n_nontripod expected 2, got %d', perexp.n_nontripod);
fprintf('Test 4 PASSED: computePerExpTCS concat across flies\n');

%% --- Test 5: computePerExpTCS with no tripod events anywhere ---------------
% All walks have candidate steps but none qualify as tripod events.
perwalk_empty = make_walkstruct({[], []}, {[], []}, [3 5], [3 5]);
perexp_empty = computePerExpTCS(perwalk_empty, struct('TCS', {[]}));
assert(isempty(perexp_empty.data), 'Test 5: data should be empty');
assert(isnan(perexp_empty.mean), 'Test 5: mean should be NaN');
assert(perexp_empty.n == 0, 'Test 5: n should be 0');
assert(perexp_empty.n_steps == 8, 'Test 5: n_steps expected 8, got %d', perexp_empty.n_steps);
assert(perexp_empty.n_nontripod == 8, 'Test 5: n_nontripod expected 8, got %d', perexp_empty.n_nontripod);
fprintf('Test 5 PASSED: experiment with no tripod events -> NaN mean, counts preserved\n');

fprintf('\n=== ALL test_TCS_aggregation TESTS PASSED ===\n');


%% ---- helpers --------------------------------------------------------------
function wfs = make_walkstruct(data_cells, speed_cells, n_steps, n_nontripod)
% Build a walkfeaturestruct-like array with only the .TCS.both fields that
% the aggregation functions read. nwalks=0 yields a valid 0-element array.
nwalks = numel(data_cells);
wfs = struct('TCS', cell(1, nwalks));
for w = 1:nwalks
    both = struct;
    both.data         = data_cells{w};
    both.event_speeds = speed_cells{w};
    both.mean         = mean(data_cells{w}, 'omitnan');
    both.std          = std(data_cells{w}, 'omitnan');
    both.n            = numel(data_cells{w});
    both.n_steps      = n_steps(w);
    both.n_nontripod  = n_nontripod(w);
    wfs(w).TCS.both = both;
end
end
