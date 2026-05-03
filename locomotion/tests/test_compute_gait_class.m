%% test_compute_gait_class.m
% Synthetic-data unit test for compute_gait_class.
%
% Coverage:
%   1. All 5 classes detected on hand-built single-frame patterns
%   2. Both tripod patterns (A and B)
%   3. All 6 tetrapod patterns (right cycle 3 + left cycle 3, including
%      the previously-missing (L1,R3) and (L3,R1))
%   4. "Other" for non-canonical patterns
%   5. Multi-fly cell input
%   6. Empty input handled
%
% Reminder: groundcontact rows are tip_pos_body order
% (RF=1, RM=2, RH=3, LH=4, LM=5, LF=6); compute_gait_class permutes
% to Mendes (L1 L2 L3 R1 R2 R3) internally.

addpath('/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis');
modpath;

fprintf('=== test_compute_gait_class ===\n\n');

k = gait_pattern_constants();

% Helper: build a 6x1 groundcontact column in tip_pos_body order from a
% Mendes-order pattern (L1 L2 L3 R1 R2 R3). The inverse of legorder_to_mendes.
mendes_to_legorder = zeros(1, 6);
mendes_to_legorder(k.legorder_to_mendes) = 1:6;
make_gc = @(mendes_pat) reshape(mendes_pat(mendes_to_legorder), 6, 1);

%% --- Test 1: tripod A and B ------------------------------------------------
gc = [make_gc([0 1 0 1 0 1]), make_gc([1 0 1 0 1 0])];   % 6 x 2
data = compute_gait_class({gc});
assert(isequal(data{1}, uint8([k.code.tripod, k.code.tripod])), ...
    'Test 1: both frames should be tripod, got %s', mat2str(data{1}));
fprintf('Test 1 PASSED: both tripod patterns classified as tripod\n');

%% --- Test 2: all 6 tetrapod patterns ---------------------------------------
% In Mendes (L1 L2 L3 R1 R2 R3) order
tet_patterns = [1 0 1 0 1 1;   % swing (L2, R1)
                1 1 0 1 0 1;   % swing (L3, R2)
                0 1 1 1 1 0;   % swing (L1, R3)  <- previously missing
                0 1 1 1 0 1;   % swing (L1, R2)
                1 0 1 1 1 0;   % swing (L2, R3)
                1 1 0 0 1 1];  % swing (L3, R1)  <- previously missing
gc = zeros(6, 6);
for i = 1:6
    gc(:, i) = make_gc(tet_patterns(i, :));
end
data = compute_gait_class({gc});
assert(all(data{1} == k.code.tetrapod), ...
    'Test 2: all 6 frames should be tetrapod, got %s', mat2str(data{1}));
fprintf('Test 2 PASSED: all 6 tetrapod patterns classified (incl. (L1,R3) and (L3,R1))\n');

%% --- Test 3: grounded (all stance) and airborne (all swing) ----------------
gc = [ones(6, 1), zeros(6, 1)];
data = compute_gait_class({gc});
assert(data{1}(1) == k.code.grounded, 'Test 3: all-1s should be grounded, got %d', data{1}(1));
assert(data{1}(2) == k.code.airborne, 'Test 3: all-0s should be airborne, got %d', data{1}(2));
fprintf('Test 3 PASSED: grounded and airborne detected\n');

%% --- Test 4: other (1 leg in swing, 5 in stance) ---------------------------
gc = make_gc([0 1 1 1 1 1]);   % only L1 in swing
data = compute_gait_class({gc});
assert(data{1} == k.code.other, 'Test 4: 5-stance pattern should be other, got %d', data{1});
fprintf('Test 4 PASSED: non-canonical pattern -> other\n');

%% --- Test 5: mixed sequence ------------------------------------------------
gc = [make_gc([0 1 0 1 0 1]), ...   % tripod A
      make_gc([1 0 1 0 1 1]), ...   % tetrapod swing(L2,R1)
      ones(6, 1), ...                 % grounded
      zeros(6, 1), ...                % airborne
      make_gc([0 1 1 1 1 1])];       % other
data = compute_gait_class({gc});
expected = uint8([k.code.tripod, k.code.tetrapod, k.code.grounded, k.code.airborne, k.code.other]);
assert(isequal(data{1}, expected), ...
    'Test 5: expected %s, got %s', mat2str(expected), mat2str(data{1}));
fprintf('Test 5 PASSED: mixed sequence classified correctly\n');

%% --- Test 6: multi-fly input -----------------------------------------------
gc1 = make_gc([0 1 0 1 0 1]);
gc2 = ones(6, 1);
data = compute_gait_class({gc1, gc2});
assert(numel(data) == 2, 'Test 6: expected 2 cells, got %d', numel(data));
assert(data{1} == k.code.tripod && data{2} == k.code.grounded, ...
    'Test 6: per-fly classification mismatch');
fprintf('Test 6 PASSED: multi-fly input handled\n');

%% --- Test 7: empty fly -----------------------------------------------------
data = compute_gait_class({zeros(6, 0)});
assert(isempty(data{1}), 'Test 7: empty input should yield empty output');
assert(isa(data{1}, 'uint8'), 'Test 7: output type should be uint8');
fprintf('Test 7 PASSED: empty input -> empty uint8\n');

%% --- Test 8: cross-check vs Alice''s 4-pattern classifier on a tetrapod -----
% The two newly-added patterns must NOT be classified as "other" anymore.
gc = [make_gc([0 1 1 1 1 0]), make_gc([1 1 0 0 1 1])];
data = compute_gait_class({gc});
assert(all(data{1} == k.code.tetrapod), ...
    'Test 8: previously-missing tetrapod patterns must classify as tetrapod, got %s', ...
    mat2str(data{1}));
fprintf('Test 8 PASSED: previously-missing tetrapod patterns no longer "other"\n');

fprintf('\nAll compute_gait_class tests PASSED.\n');
