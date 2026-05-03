%% test_compute_TCS.m
% Synthetic-data unit test for compute_TCS (revised algorithm).
%
% Coverage:
%   1. Perfect tripod step (3 legs swing simultaneously) -> TCS = 1
%   2. Sequential swings, no overlap -> NOT a tripod event (n=0, n_nontripod=1)
%   3. Partial overlap with t2 >= 1 -> valid event with TCS = ts/ti
%   4. Multi-event walk: mix of tripod and non-tripod steps; tripod fraction
%   5. Both-tripod aggregation
%   6. Walk window respect
%   7. Restriction C (Hilbert-valid window) excludes events outside intersection
%   8. Per-event speed = mean(velmag_ctr) over t1 envelope
%   9. step_idx records the lead-leg's step index for each event
%
% Tripod groupings (compute_TCS.m):
%   tripod_LM = {LM=5 (lead), RF=1, RH=3}
%   tripod_RM = {RM=2 (lead), LF=6, LH=4}
%
% All bout indices follow the detect_bouts convention:
%   start_indices: inclusive, end_indices: exclusive (first frame after).

addpath('/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis');
modpath;

fprintf('=== test_compute_TCS ===\n\n');

%% --- Test 1: perfect tripod, single LM step --------------------------------
% Walk window 1..51. LM has one full step cycle (stance [1,10), swing [10,20),
% next stance starts at 20). RF and RH swing simultaneously [10,20).
% Entire walk has the same Hilbert peaks for all 3 legs (frames 5 and 15).
walk_t0 = 1; walk_t1 = 51;

bd = make_empty_bd();
% LM (lead): stance ends at 10, swing [10,20), next stance at 20
bd = add_step(bd, 5, 1, 10, 20);     % step 1: stance [1,10), swing [10,20), end 20
% RF: matching swing [10, 20)
bd.perlimb(1).swing.start_indices(end+1,1) = 10;
bd.perlimb(1).swing.end_indices(end+1,1)   = 20;
% RH: matching swing [10, 20)
bd.perlimb(3).swing.start_indices(end+1,1) = 10;
bd.perlimb(3).swing.end_indices(end+1,1)   = 20;

% Wide Hilbert windows for all 6 legs (peaks at [5, 15, 25, 35, 45])
loctall = make_peak_cells(6, [5 15 25 35 45]);
locball = make_peak_cells(6, [5 15 25 35 45]);

speed_data = ones(1, 100) * 10;  % constant 10 mm/s

TCS = compute_TCS(bd, walk_t0, walk_t1, loctall, locball, speed_data);

assert(TCS.tripod_LM.n == 1, 'Test 1: expected n=1 valid LM event, got %d', TCS.tripod_LM.n);
assert(TCS.tripod_LM.n_nontripod == 0, 'Test 1: expected n_nontripod=0, got %d', TCS.tripod_LM.n_nontripod);
assert(abs(TCS.tripod_LM.data(1) - 1) < 1e-12, 'Test 1: TCS expected 1, got %g', TCS.tripod_LM.data(1));
assert(TCS.tripod_LM.event_times(1) == 1, 'Test 1: event_time should be step start = 1, got %d', TCS.tripod_LM.event_times(1));
assert(TCS.tripod_LM.step_idx(1) == 1, 'Test 1: step_idx should be 1, got %d', TCS.tripod_LM.step_idx(1));
assert(abs(TCS.tripod_LM.event_speeds(1) - 10) < 1e-12, 'Test 1: speed expected 10, got %g', TCS.tripod_LM.event_speeds(1));
fprintf('Test 1 PASSED: perfect tripod LM step -> TCS=1, time=1, step_idx=1, speed=10\n');

%% --- Test 2: sequential swings, no overlap ---------------------------------
% LM swings [10,15), RF [16,21), RH [22,27) -- no all-3-overlap -> not a tripod step
walk_t0 = 1; walk_t1 = 51;
bd = make_empty_bd();
bd = add_step(bd, 5, 1, 10, 28);     % LM step covering [1,28); swing [10,15)
bd.perlimb(5).swing.end_indices(end) = 15;   % override to make swing [10,15)
bd.perlimb(1).swing.start_indices(end+1,1) = 16;
bd.perlimb(1).swing.end_indices(end+1,1)   = 21;
bd.perlimb(3).swing.start_indices(end+1,1) = 22;
bd.perlimb(3).swing.end_indices(end+1,1)   = 27;

loctall = make_peak_cells(6, [5 15 25 35 45]);
locball = make_peak_cells(6, [5 15 25 35 45]);
speed_data = ones(1, 100) * 10;

TCS = compute_TCS(bd, 1, 51, loctall, locball, speed_data);
assert(TCS.tripod_LM.n == 0, 'Test 2: expected n=0 (no valid tripod events), got %d', TCS.tripod_LM.n);
assert(TCS.tripod_LM.n_steps == 1, 'Test 2: expected 1 candidate step, got %d', TCS.tripod_LM.n_steps);
assert(TCS.tripod_LM.n_nontripod == 1, 'Test 2: expected n_nontripod=1, got %d', TCS.tripod_LM.n_nontripod);
fprintf('Test 2 PASSED: sequential swings -> not tripod (n_steps=1, n=0, n_nontripod=1)\n');

%% --- Test 3: partial overlap with t2 >= 1 ----------------------------------
% LM swing [10,20), RF [12,18), RH [14,22) -- all-3 overlap [14,18) -> ts=4, ti=12 -> TCS=1/3
walk_t0 = 1; walk_t1 = 51;
bd = make_empty_bd();
bd = add_step(bd, 5, 1, 10, 23);
bd.perlimb(5).swing.end_indices(end) = 20;
bd.perlimb(1).swing.start_indices(end+1,1) = 12;
bd.perlimb(1).swing.end_indices(end+1,1)   = 18;
bd.perlimb(3).swing.start_indices(end+1,1) = 14;
bd.perlimb(3).swing.end_indices(end+1,1)   = 22;

loctall = make_peak_cells(6, [5 15 25 35 45]);
locball = make_peak_cells(6, [5 15 25 35 45]);
speed_data = ones(1, 100) * 10;

TCS = compute_TCS(bd, 1, 51, loctall, locball, speed_data);
assert(TCS.tripod_LM.n == 1, 'Test 3: expected n=1, got %d', TCS.tripod_LM.n);
expected_tcs = 4 / 12;
assert(abs(TCS.tripod_LM.data(1) - expected_tcs) < 1e-12, ...
    'Test 3: TCS expected %g, got %g', expected_tcs, TCS.tripod_LM.data(1));
fprintf('Test 3 PASSED: partial overlap -> TCS = 4/12 = %.4f\n', expected_tcs);

%% --- Test 4: multi-event walk: 3 LM steps, mixed valid/invalid -------------
% Step 1: perfect overlap [10,20) -> TCS=1
% Step 2: partial overlap [40,50)/[40,50)/[40,60) -> TCS = 10/20 = 0.5
% Step 3: sequential -> not a tripod event
walk_t0 = 1; walk_t1 = 100;
bd = make_empty_bd();
% LM steps: [1,30), [30,70), [70,100)
bd.perlimb(5).step.start_indices = [1; 30; 70];
bd.perlimb(5).step.end_indices   = [30; 70; 100];
bd.perlimb(5).swing.start_indices = [10; 40; 80];
bd.perlimb(5).swing.end_indices   = [20; 50; 85];
% RF
bd.perlimb(1).swing.start_indices = [10; 40; 86];
bd.perlimb(1).swing.end_indices   = [20; 50; 91];
% RH
bd.perlimb(3).swing.start_indices = [10; 40; 92];
bd.perlimb(3).swing.end_indices   = [20; 60; 97];

loctall = make_peak_cells(6, [5 25 50 75 95]);
locball = make_peak_cells(6, [5 25 50 75 95]);
speed_data = ones(1, 100) * 10;

TCS = compute_TCS(bd, 1, 100, loctall, locball, speed_data);
assert(TCS.tripod_LM.n_steps == 3, 'Test 4: expected 3 candidate steps, got %d', TCS.tripod_LM.n_steps);
assert(TCS.tripod_LM.n == 2, 'Test 4: expected n=2 valid tripod events, got %d', TCS.tripod_LM.n);
assert(TCS.tripod_LM.n_nontripod == 1, 'Test 4: expected n_nontripod=1, got %d', TCS.tripod_LM.n_nontripod);
expected_data = [1, 0.5];
assert(all(abs(TCS.tripod_LM.data - expected_data) < 1e-12), ...
    'Test 4: data expected %s, got %s', mat2str(expected_data), mat2str(TCS.tripod_LM.data));
expected_step_idx = [1, 2];
assert(isequal(TCS.tripod_LM.step_idx, expected_step_idx), ...
    'Test 4: step_idx expected %s, got %s', mat2str(expected_step_idx), mat2str(TCS.tripod_LM.step_idx));
fprintf('Test 4 PASSED: 3 steps, 2 tripod (TCS=[1,0.5]), 1 non-tripod, step_idx=[1,2]\n');

%% --- Test 5: both-tripod aggregation ---------------------------------------
% LM tripod has 1 perfect event, RM tripod has 1 perfect event. Both n=1, mean=1.
walk_t0 = 1; walk_t1 = 51;
bd = make_empty_bd();
% LM perfect step
bd = add_step(bd, 5, 1, 10, 25);     % LM step [1,25), swing [10,25) per add_step
bd.perlimb(5).swing.end_indices(end) = 20;     % override to make LM swing [10,20)
bd.perlimb(1).swing.start_indices(end+1,1) = 10; bd.perlimb(1).swing.end_indices(end+1,1) = 20;  % RF
bd.perlimb(3).swing.start_indices(end+1,1) = 10; bd.perlimb(3).swing.end_indices(end+1,1) = 20;  % RH
% RM perfect step
bd = add_step(bd, 2, 25, 35, 50);    % RM step [25,50), swing [35,50) per add_step
bd.perlimb(2).swing.end_indices(end) = 45;     % override to make RM swing [35,45)
bd.perlimb(6).swing.start_indices(end+1,1) = 35; bd.perlimb(6).swing.end_indices(end+1,1) = 45;  % LF
bd.perlimb(4).swing.start_indices(end+1,1) = 35; bd.perlimb(4).swing.end_indices(end+1,1) = 45;  % LH

loctall = make_peak_cells(6, [5 20 30 45]);
locball = make_peak_cells(6, [5 20 30 45]);
speed_data = ones(1, 100) * 10;

TCS = compute_TCS(bd, 1, 51, loctall, locball, speed_data);
assert(TCS.both.n == 2, 'Test 5: expected both.n=2, got %d', TCS.both.n);
assert(abs(TCS.both.mean - 1) < 1e-12, 'Test 5: expected both.mean=1, got %g', TCS.both.mean);
assert(isequal(TCS.both.tripod_id, [1 2]), 'Test 5: expected tripod_id=[1 2], got %s', mat2str(TCS.both.tripod_id));
fprintf('Test 5 PASSED: both-tripod aggregation, n=2, mean=1.0, tripod_id=[1 2]\n');

%% --- Test 6: walk window respect -------------------------------------------
% LM steps at trajectory frames [5,30), [30,55), [80,100). Walk window [10, 50).
% Only step starting at 30 has start in walk window.
walk_t0 = 10; walk_t1 = 50;
bd = make_empty_bd();
bd.perlimb(5).step.start_indices = [5; 30; 80];
bd.perlimb(5).step.end_indices   = [30; 55; 100];
bd.perlimb(5).swing.start_indices = [20; 40; 90];
bd.perlimb(5).swing.end_indices   = [25; 45; 95];
bd.perlimb(1).swing.start_indices = [20; 40; 90];
bd.perlimb(1).swing.end_indices   = [25; 45; 95];
bd.perlimb(3).swing.start_indices = [20; 40; 90];
bd.perlimb(3).swing.end_indices   = [25; 45; 95];

loctall = make_peak_cells(6, [3 23 43 70 95]);
locball = make_peak_cells(6, [3 23 43 70 95]);
speed_data = ones(1, 100) * 10;

TCS = compute_TCS(bd, walk_t0, walk_t1, loctall, locball, speed_data);
assert(TCS.tripod_LM.n_steps == 1, 'Test 6: expected 1 candidate step, got %d', TCS.tripod_LM.n_steps);
assert(TCS.tripod_LM.n == 1, 'Test 6: expected n=1, got %d', TCS.tripod_LM.n);
assert(TCS.tripod_LM.step_idx(1) == 2, 'Test 6: step_idx expected 2, got %d', TCS.tripod_LM.step_idx(1));
fprintf('Test 6 PASSED: walk window filters steps (only step 2 of 3 in window)\n');

%% --- Test 7: Restriction C excludes events outside Hilbert-valid window ----
% LM steps at frames [1,30), [30,60), [60,100). All 3 perfect tripod events.
% LM peaks at [25, 55] (so valid window = [25, 55]). RF and RH have wide windows.
% Step 1 (start=1) outside; step 2 (start=30) inside; step 3 (start=60) outside.
walk_t0 = 1; walk_t1 = 100;
bd = make_empty_bd();
bd.perlimb(5).step.start_indices = [1; 30; 60];
bd.perlimb(5).step.end_indices   = [30; 60; 100];
bd.perlimb(5).swing.start_indices = [10; 40; 70];
bd.perlimb(5).swing.end_indices   = [20; 50; 80];
bd.perlimb(1).swing.start_indices = [10; 40; 70];
bd.perlimb(1).swing.end_indices   = [20; 50; 80];
bd.perlimb(3).swing.start_indices = [10; 40; 70];
bd.perlimb(3).swing.end_indices   = [20; 50; 80];

% Only LM has a narrow Hilbert window
loctall = make_peak_cells(6, [3 50 95]);
locball = make_peak_cells(6, [3 50 95]);
loctall{5} = [0 0; 25 55];   % LM peaks at local frames 25 and 55
locball{5} = [0 0; 25 55];

speed_data = ones(1, 100) * 10;

TCS = compute_TCS(bd, 1, 100, loctall, locball, speed_data);
assert(TCS.tripod_LM.n_steps == 1, 'Test 7: expected only 1 candidate step inside Hilbert window, got %d', TCS.tripod_LM.n_steps);
assert(TCS.tripod_LM.n == 1, 'Test 7: expected 1 valid tripod event, got %d', TCS.tripod_LM.n);
assert(TCS.tripod_LM.step_idx(1) == 2, 'Test 7: only step 2 should be in Hilbert window, got step_idx=%d', TCS.tripod_LM.step_idx(1));
fprintf('Test 7 PASSED: Restriction C excludes steps outside Hilbert-valid window\n');

%% --- Test 8: per-event speed averaged over t1 envelope ---------------------
% LM swing [10,20), RF [12,18), RH [14,22).
% t1 envelope: [min(10,12,14), max(20,18,22)) = [10, 22) -> indices 10..21 inclusive
% speed_data ramps 1..100. mean(10..21) = (10+21)/2 = 15.5
walk_t0 = 1; walk_t1 = 51;
bd = make_empty_bd();
bd = add_step(bd, 5, 1, 10, 25);
bd.perlimb(5).swing.end_indices(end) = 20;
bd.perlimb(1).swing.start_indices(end+1,1) = 12; bd.perlimb(1).swing.end_indices(end+1,1) = 18;
bd.perlimb(3).swing.start_indices(end+1,1) = 14; bd.perlimb(3).swing.end_indices(end+1,1) = 22;

loctall = make_peak_cells(6, [5 25 45]);
locball = make_peak_cells(6, [5 25 45]);
speed_data = 1:100;

TCS = compute_TCS(bd, 1, 51, loctall, locball, speed_data);
expected_speed = mean(10:21);
assert(abs(TCS.tripod_LM.event_speeds(1) - expected_speed) < 1e-12, ...
    'Test 8: speed expected %g (mean of frames 10..21), got %g', ...
    expected_speed, TCS.tripod_LM.event_speeds(1));
fprintf('Test 8 PASSED: speed averaged over t1 envelope = mean(10:21) = %.1f\n', expected_speed);

fprintf('\n=== ALL TESTS PASSED ===\n');


%% ---- helpers at end of script (MATLAB requires this order) ----

function bd = make_empty_bd()
% Empty per-fly bout struct with 6 legs and swing/stance/step subfields.
bd = struct();
for L = 1:6
    bd.perlimb(L).swing.start_indices  = zeros(0,1);
    bd.perlimb(L).swing.end_indices    = zeros(0,1);
    bd.perlimb(L).stance.start_indices = zeros(0,1);
    bd.perlimb(L).stance.end_indices   = zeros(0,1);
    bd.perlimb(L).step.start_indices   = zeros(0,1);
    bd.perlimb(L).step.end_indices     = zeros(0,1);
end
end

function bd = add_step(bd, leg, step_start, swing_start, step_end)
% Append one step to leg `leg`: stance [step_start, swing_start),
% swing [swing_start, step_end), step [step_start, step_end).
bd.perlimb(leg).step.start_indices(end+1,1)   = step_start;
bd.perlimb(leg).step.end_indices(end+1,1)     = step_end;
bd.perlimb(leg).swing.start_indices(end+1,1)  = swing_start;
bd.perlimb(leg).swing.end_indices(end+1,1)    = step_end;
end

function pk = make_peak_cells(nlegs, frame_list)
% Build a {1 x nlegs} cell where each entry has the same set of peak frames
% (used for Hilbert-window restriction tests; values aren't used).
pk = cell(1, nlegs);
for L = 1:nlegs
    pk{L} = [zeros(1, numel(frame_list)); frame_list(:)'];
end
end
