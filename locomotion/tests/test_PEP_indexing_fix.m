%% test_PEP_indexing_fix.m
% Synthetic-data test of two related claims:
%   1. detect_bouts() returns INCLUSIVE start_indices and EXCLUSIVE end_indices
%      (end_indices points to the FIRST FRAME AFTER the bout — i.e., the first
%      non-stance / first swing frame).
%   2. computeStepFeatures.m indexes AEP/PEP per the velocity-interval convention:
%        velocity[t] = |pos[t+1]-pos[t]| (forward difference), so gc[t]=1 means the
%        foot is planted over interval [t, t+1].
%        - AEP comes from frame stance_t0s  (touchdown, first in-contact frame)
%        - PEP comes from frame stance_t1s  (lift-off: foot is stationary through
%          stance_t1s and first moves over [stance_t1s, stance_t1s+1])
%
% History: PEP was previously pos[stance_t1s-1], which is one frame too early under
% the velocity-interval convention (the foot is still planted at stance_t1s). The
% even-older code used pos[stance_t1s] for a different (incorrect) reason; this test
% pins the CURRENT, correct definition PEP = pos[stance_t1s].

addpath('/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis');
modpath;

fprintf('=== test_PEP_indexing_fix ===\n\n');

%% --- Test 1: detect_bouts convention ----------------------------------------
% Synthetic groundcontact: bouts (1s) at frames 3-7 and 12-15.
% Frames:    1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
% Stance(1)= 0 0 1 1 1 1 1 0 0  0  0  1  1  1  1  0
gc = [0 0 1 1 1 1 1 0 0 0 0 1 1 1 1 0]';

[start_indices, end_indices] = detect_bouts(gc);

% Expected: bout 1 starts at frame 3, ends (exclusive) at frame 8 (first 0 after).
%           bout 2 starts at frame 12, ends (exclusive) at frame 16.
expected_starts = [3; 12];
expected_ends   = [8; 16];

assert(isequal(start_indices(:), expected_starts), ...
    'detect_bouts start_indices wrong: got %s, expected %s', ...
    mat2str(start_indices(:)'), mat2str(expected_starts'));
assert(isequal(end_indices(:), expected_ends), ...
    'detect_bouts end_indices wrong: got %s, expected %s', ...
    mat2str(end_indices(:)'), mat2str(expected_ends'));

% Confirm semantics: gc(start) == 1 (first stance frame)
%                    gc(end-1) == 1 (last stance frame)
%                    gc(end)   == 0 (first swing frame after stance)
for i = 1:numel(start_indices)
    assert(gc(start_indices(i)) == 1, 'gc(start) should be 1 (in-stance)');
    assert(gc(end_indices(i) - 1) == 1, 'gc(end-1) should be 1 (last stance frame)');
    assert(gc(end_indices(i)) == 0, 'gc(end) should be 0 (first swing frame, EXCLUSIVE)');
end
fprintf('Test 1 PASSED: detect_bouts returns inclusive start, exclusive end.\n');
fprintf('  Bout 1: frames 3-7 stance, end_index = 8 (first swing frame)\n');
fprintf('  Bout 2: frames 12-15 stance, end_index = 16 (first swing frame)\n\n');

%% --- Test 2: computeStepFeatures AEP/PEP indexing ---------------------------
% Build a synthetic tip_pos_body where each frame's leg position is identifiable:
%   tip_pos_body(limb, 1, t) = t        (x = frame number)
%   tip_pos_body(limb, 2, t) = t * 10   (y = frame * 10)
% This way the value of PEP/AEP tells us exactly which frame was indexed.

nlimbs   = 6;
nframes  = 30;
nlandmrk = 6;
limb     = 1;

tip_pos_body = nan(nlimbs, 2, nframes);
for t = 1:nframes
    for L = 1:nlimbs
        tip_pos_body(L, 1, t) = t;
        tip_pos_body(L, 2, t) = t * 10;
    end
end

% Define stance/step bouts using the detect_bouts convention from Test 1
stance_t0s = [3; 12];
stance_t1s = [8; 16];           % EXCLUSIVE
step_t0s   = [3; 12];            % step starts at touchdown
step_t1s   = [12; 16];           % step ends at next touchdown (exclusive)

% Minimal mock inputs so the function runs end-to-end
fly = 1;
trx = struct();
trx(fly).a = ones(1, nframes) * 0.5;        % body half-length; mean*4 = 2 px
trx(fly).off = 0;

% APT data needs to exist for the lab-frame `length_px` computation
% (lines 90-100 of computeStepFeatures); shape: legtip_landmarknums rows
aptdata = struct();
aptdata.pTrk = cell(1, fly);
aptdata.pTrk{fly} = nan(nlandmrk, 2, nframes);
for t = 1:nframes
    aptdata.pTrk{fly}(:, 1, t) = t;
    aptdata.pTrk{fly}(:, 2, t) = t * 10;
end
legtip_landmarknums = 1:nlandmrk;

currfly_timestamps = (0:nframes-1)' / 150;   % 150 fps

%% Call the function
boutfeatures = computeStepFeatures(fly, trx, aptdata, tip_pos_body, ...
    legtip_landmarknums, limb, step_t0s, step_t1s, stance_t0s, stance_t1s, ...
    currfly_timestamps);

%% Check AEP — should be at stance_t0s frames (3 and 12)
expected_AEP = [3, 12; 30, 120];     % [x; y] for each bout
assert(isequal(boutfeatures.AEP, expected_AEP), ...
    'AEP wrong: got %s, expected %s', ...
    mat2str(boutfeatures.AEP), mat2str(expected_AEP));
fprintf('Test 2a PASSED: AEP at stance_t0s = [3, 12] (touchdown frames).\n');
fprintf('  AEP = %s (x), %s (y)\n', ...
    mat2str(boutfeatures.AEP(1,:)), mat2str(boutfeatures.AEP(2,:)));

%% Check PEP — should be at stance_t1s = [8, 16] (lift-off, first swing index)
%   prior (off-by-one) definition was stance_t1s-1 = [7, 15]
expected_PEP          = [8, 16; 80, 160];
expected_PEP_oldoff1  = [7, 15; 70, 150];

if isequal(boutfeatures.PEP, expected_PEP)
    fprintf('Test 2b PASSED: PEP at stance_t1s = [8, 16] (lift-off position).\n');
    fprintf('  PEP = %s (x), %s (y)\n', ...
        mat2str(boutfeatures.PEP(1,:)), mat2str(boutfeatures.PEP(2,:)));
elseif isequal(boutfeatures.PEP, expected_PEP_oldoff1)
    error(['Test 2b FAILED: PEP indexed at stance_t1s-1 = [7, 15] (one frame too ', ...
           'early). The PEP fix in computeStepFeatures.m was NOT applied.']);
else
    error('Test 2b FAILED: PEP = %s; expected %s (correct) or %s (old off-by-one).', ...
        mat2str(boutfeatures.PEP), mat2str(expected_PEP), ...
        mat2str(expected_PEP_oldoff1));
end

%% Check amplitude — ‖PEP - AEP‖ in pixels, then in body lengths
% bout 1 amplitude = ‖[8-3; 80-30]‖ = ‖[5; 50]‖ = sqrt(25+2500) ≈ 50.25
% bout 2 amplitude = ‖[16-12; 160-120]‖ = ‖[4; 40]‖ = sqrt(16+1600) ≈ 40.20
expected_amplitude_px = [hypot(5, 50); hypot(4, 40)];
assert(all(abs(boutfeatures.amplitude_px(:) - expected_amplitude_px) < 1e-9), ...
    'amplitude_px wrong: got %s, expected %s', ...
    mat2str(boutfeatures.amplitude_px(:)'), mat2str(expected_amplitude_px'));
fprintf('Test 2c PASSED: amplitude_px matches expected (PEP = stance_t1s).\n');
fprintf('  amplitude_px = %s\n', mat2str(boutfeatures.amplitude_px));

%% Check step_direction — atan2(PEP_y - AEP_y, PEP_x - AEP_x)
% bout 1 direction = atan2(80-30, 8-3)   = atan2(50, 5)
% bout 2 direction = atan2(160-120, 16-12) = atan2(40, 4)
expected_step_dir = [atan2(50, 5), atan2(40, 4)];
assert(all(abs(boutfeatures.step_direction(:) - expected_step_dir(:)) < 1e-9), ...
    'step_direction wrong: got %s, expected %s', ...
    mat2str(boutfeatures.step_direction), mat2str(expected_step_dir));
fprintf('Test 2d PASSED: step_direction matches expected (PEP = stance_t1s).\n\n');

fprintf('=== ALL TESTS PASSED ===\n');
