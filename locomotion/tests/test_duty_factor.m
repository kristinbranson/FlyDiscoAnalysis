%% test_duty_factor.m
% Unit test for duty_factor computation in computeStepFeatures.
%
% Coverage:
%   1. Known values: stance=60ms, swing=40ms -> duty_factor = 0.6
%   2. Edge case: step duration = 0 -> duty_factor = NaN
%   3. Edge case: all stance (swing=0) -> duty_factor = 1.0
%   4. Edge case: all swing (stance=0) -> duty_factor = 0.0
%   5. Multiple steps: verify vectorized computation
%   6. Verify duty_factor propagates through combineStepMetrics

addpath('/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis');
modpath;

fprintf('=== test_duty_factor ===\n\n');
npass = 0;
nfail = 0;

%% --- Test 1: Known values (stance=60ms, swing=40ms) -----------------------
fprintf('Test 1: Known values (stance=60ms, swing=40ms, duty_factor=0.6)...\n');
% fps=1000 so 1 frame = 1ms. Step = [1,101), stance = [1,61), swing = [61,101).
timestamps = (0:200) * 0.001; % seconds, so *1000 = ms
step_t0s = 1;
step_t1s = 101;  % exclusive: step is frames 1..100 = 100 frames = 100ms
stance_t0s = 1;
stance_t1s = 61;  % exclusive: stance is frames 1..60 = 60 frames = 60ms

[~, dur_time] = computeBoutDurations(step_t0s, step_t1s, timestamps);
[~, stance_dur_time] = computeBoutDurations(stance_t0s, stance_t1s, timestamps);
duty_factor = stance_dur_time ./ dur_time;
duty_factor(dur_time == 0) = NaN;

expected = 0.6;
if abs(duty_factor - expected) < 1e-10
    fprintf('  PASS: duty_factor = %.4f (expected %.4f)\n', duty_factor, expected);
    npass = npass + 1;
else
    fprintf('  FAIL: duty_factor = %.4f (expected %.4f)\n', duty_factor, expected);
    nfail = nfail + 1;
end

%% --- Test 2: Zero step duration -> NaN ------------------------------------
fprintf('Test 2: Zero step duration -> NaN...\n');
step_t0s = 50;
step_t1s = 50;  % same frame -> 0 duration
stance_t0s = 50;
stance_t1s = 50;

[~, dur_time] = computeBoutDurations(step_t0s, step_t1s, timestamps);
[~, stance_dur_time] = computeBoutDurations(stance_t0s, stance_t1s, timestamps);
duty_factor = stance_dur_time ./ dur_time;
duty_factor(dur_time == 0) = NaN;

if isnan(duty_factor)
    fprintf('  PASS: duty_factor = NaN\n');
    npass = npass + 1;
else
    fprintf('  FAIL: duty_factor = %.4f (expected NaN)\n', duty_factor);
    nfail = nfail + 1;
end

%% --- Test 3: All stance (swing=0) -> duty_factor = 1.0 -------------------
fprintf('Test 3: All stance (swing=0) -> duty_factor = 1.0...\n');
step_t0s = 1;
step_t1s = 51;   % 50 frames
stance_t0s = 1;
stance_t1s = 51;  % stance fills entire step

[~, dur_time] = computeBoutDurations(step_t0s, step_t1s, timestamps);
[~, stance_dur_time] = computeBoutDurations(stance_t0s, stance_t1s, timestamps);
duty_factor = stance_dur_time ./ dur_time;
duty_factor(dur_time == 0) = NaN;

if abs(duty_factor - 1.0) < 1e-10
    fprintf('  PASS: duty_factor = %.4f\n', duty_factor);
    npass = npass + 1;
else
    fprintf('  FAIL: duty_factor = %.4f (expected 1.0)\n', duty_factor);
    nfail = nfail + 1;
end

%% --- Test 4: All swing (stance=0) -> duty_factor = 0.0 -------------------
fprintf('Test 4: All swing (stance=0) -> duty_factor = 0.0...\n');
step_t0s = 1;
step_t1s = 51;   % 50 frames
stance_t0s = 1;
stance_t1s = 1;   % 0 frames of stance

[~, dur_time] = computeBoutDurations(step_t0s, step_t1s, timestamps);
[~, stance_dur_time] = computeBoutDurations(stance_t0s, stance_t1s, timestamps);
duty_factor = stance_dur_time ./ dur_time;
duty_factor(dur_time == 0) = NaN;

if abs(duty_factor - 0.0) < 1e-10
    fprintf('  PASS: duty_factor = %.4f\n', duty_factor);
    npass = npass + 1;
else
    fprintf('  FAIL: duty_factor = %.4f (expected 0.0)\n', duty_factor);
    nfail = nfail + 1;
end

%% --- Test 5: Multiple steps, vectorized ----------------------------------
fprintf('Test 5: Multiple steps (3 steps with different duty factors)...\n');
% Step 1: frames 1-100 (100ms), stance 1-60 (60ms) -> DF=0.6
% Step 2: frames 101-200 (100ms), stance 101-180 (80ms) -> DF=0.8
% Step 3: frames 201-251 (50ms), stance 201-226 (25ms) -> DF=0.5
step_t0s = [1, 101, 201];
step_t1s = [101, 201, 251];
stance_t0s = [1, 101, 201];
stance_t1s = [61, 181, 226];

timestamps_long = (0:300) * 0.001;
[~, dur_time] = computeBoutDurations(step_t0s, step_t1s, timestamps_long);
[~, stance_dur_time] = computeBoutDurations(stance_t0s, stance_t1s, timestamps_long);
duty_factor = stance_dur_time ./ dur_time;
duty_factor(dur_time == 0) = NaN;

expected = [0.6, 0.8, 0.5];
if all(abs(duty_factor - expected) < 1e-10)
    fprintf('  PASS: duty_factor = [%.4f, %.4f, %.4f]\n', duty_factor);
    npass = npass + 1;
else
    fprintf('  FAIL: duty_factor = [%.4f, %.4f, %.4f] (expected [0.6, 0.8, 0.5])\n', duty_factor);
    nfail = nfail + 1;
end

%% --- Test 6: Aggregation through combineStepMetrics -----------------------
fprintf('Test 6: Verify duty_factor in combineStepMetrics stepfeatures list...\n');
% Check that 'duty_factor' is in the stepfeatures list used by combineStepMetrics.
% We do this by reading the source and checking the list.
src = fileread(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'LimbBoutAnalyzer.m'));
if contains(src, '''duty_factor''') && contains(src, 'stepfeatures = {')
    % Extract the stepfeatures cell array definition
    idx = strfind(src, 'stepfeatures = {');
    % Check duty_factor appears in the stepfeatures definition
    chunk = src(idx(1):min(idx(1)+500, end));
    if contains(chunk, 'duty_factor')
        fprintf('  PASS: duty_factor found in combineStepMetrics stepfeatures list\n');
        npass = npass + 1;
    else
        fprintf('  FAIL: duty_factor not in stepfeatures list\n');
        nfail = nfail + 1;
    end
else
    fprintf('  FAIL: could not find stepfeatures list in LimbBoutAnalyzer.m\n');
    nfail = nfail + 1;
end

%% --- Test 7: Verify duty_factor in step_geom_scalar -----------------------
fprintf('Test 7: Verify duty_factor in buildWalkStructForCondition step_geom_scalar...\n');
if contains(src, 'step_geom_scalar')
    idx = strfind(src, 'step_geom_scalar = {');
    chunk = src(idx(1):min(idx(1)+500, end));
    if contains(chunk, 'duty_factor')
        fprintf('  PASS: duty_factor found in step_geom_scalar list\n');
        npass = npass + 1;
    else
        fprintf('  FAIL: duty_factor not in step_geom_scalar list\n');
        nfail = nfail + 1;
    end
else
    fprintf('  FAIL: could not find step_geom_scalar in LimbBoutAnalyzer.m\n');
    nfail = nfail + 1;
end

%% --- Test 8: computeStepFeatures end-to-end with synthetic trx ------------
fprintf('Test 8: computeStepFeatures end-to-end with synthetic data...\n');
% Build minimal trx and aptdata structs to call computeStepFeatures directly.
nframes = 300;
timestamps_syn = (0:nframes-1) * 0.001; % 1kHz

% Minimal trx struct
trx_syn = struct();
trx_syn.a = 10 * ones(1, nframes); % body length = 10*4 = 40 px
trx_syn.firstframe = 1;
trx_syn.endframe = nframes;
trx_syn.off = 0;

% tip_pos_body: 6 limbs x 2 coords x nframes
tip_pos_body_syn = randn(6, 2, nframes);

% aptdata with pTrk
aptdata_syn = struct();
aptdata_syn.pTrk = {randn(18, 2, nframes)};  % 18 keypoints, 2 coords

legtip_landmarknums = 12:17;

% 2 steps:
% Step 1: [1,101), stance [1,61), swing = step_dur - stance_dur
% Step 2: [101,201), stance [101,151)
step_t0s_syn = [1, 101];
step_t1s_syn = [101, 201];
stance_t0s_syn = [1, 101];
stance_t1s_syn = [61, 151];

fly = 1;
limb = 1;
[~, stance_dur_syn] = computeBoutDurations(stance_t0s_syn, stance_t1s_syn, timestamps_syn);
boutfeatures = computeStepFeatures(fly, trx_syn, aptdata_syn, tip_pos_body_syn, ...
    legtip_landmarknums, limb, step_t0s_syn, step_t1s_syn, stance_t0s_syn, stance_t1s_syn, timestamps_syn, stance_dur_syn);

if isfield(boutfeatures, 'duty_factor')
    expected_df = [0.6, 0.5];
    if all(abs(boutfeatures.duty_factor - expected_df) < 1e-10)
        fprintf('  PASS: boutfeatures.duty_factor = [%.4f, %.4f]\n', boutfeatures.duty_factor);
        npass = npass + 1;
    else
        fprintf('  FAIL: boutfeatures.duty_factor = [%.4f, %.4f] (expected [0.6, 0.5])\n', boutfeatures.duty_factor);
        nfail = nfail + 1;
    end
else
    fprintf('  FAIL: duty_factor field not in boutfeatures\n');
    nfail = nfail + 1;
end

%% Summary
fprintf('\n=== SUMMARY: %d passed, %d failed ===\n', npass, nfail);
if nfail > 0
    error('test_duty_factor: %d test(s) failed', nfail);
end
