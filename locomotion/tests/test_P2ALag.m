function test_P2ALag
% test_P2ALag - Unit tests for computeP2ALag on synthetic bout data.
%
% Builds synthetic swing/stance bouts with known onset offsets and verifies
% the posterior-to-anterior onset lags (ms). 1 ms per frame timestamps so
% lag(ms) == frame difference, making expected values exact.
%
% Covers: known lag (swing->swing and touchdown->swing), zero lag,
% forward-window behavior on a "reversed" wave, single onset, missing match.
%
% Run: test_P2ALag   (errors on first failed assertion; prints PASS summary)

% computeP2ALag + findClosestSteps live one dir up (locomotion/)
addpath(fileparts(fileparts(mfilename('fullpath'))));

legorder = {'RF','RM','RH','LH','LM','LF'};   %#ok<NASGU> % 1..6 reference
dt = 0.001;                                    % 1 ms per frame
N  = 500;
timestamps = (0:N-1) * dt;                     % seconds
walk_t0 = 1; walk_t1 = 400;
fly = 1;
tol = 1e-9;
swingdur = 5;                                  % frames

npass = 0;

%% Test 1: known posterior-leads lag, swing->swing and touchdown->swing
% Right side RH(3)->RM(2)->RF(1), each anterior leg lifts off 10 frames after
% the posterior. Left side LH(4)->LM(5)->LF(6) offset by 5 to differ.
cyc = [0 50 100 150];
onsets = cell(1,6);
onsets{3} = 10 + cyc;  onsets{2} = 20 + cyc;  onsets{1} = 30 + cyc;   % R: H,M,F
onsets{4} = 15 + cyc;  onsets{5} = 25 + cyc;  onsets{6} = 35 + cyc;   % L: H,M,F
bd = makeboutdata(onsets, swingdur);

lift = computeP2ALag(bd, walk_t0, walk_t1, timestamps, fly, 'swing');
assert(approx(lift.RH_to_RM.mean, 10, tol), 'T1 RH_to_RM liftoff lag');
assert(approx(lift.RM_to_RF.mean, 10, tol), 'T1 RM_to_RF liftoff lag');
assert(approx(lift.LH_to_LM.mean, 10, tol), 'T1 LH_to_LM liftoff lag');
assert(approx(lift.LM_to_LF.mean, 10, tol), 'T1 LM_to_LF liftoff lag');
assert(lift.RH_to_RM.n == 3, 'T1 RH_to_RM n');
% aggregates: concatenate-then-mean
assert(approx(lift.H_to_M.mean, 10, tol) && lift.H_to_M.n == 6, 'T1 H_to_M');
assert(approx(lift.M_to_F.mean, 10, tol) && lift.M_to_F.n == 6, 'T1 M_to_F');
assert(approx(lift.all.mean, 10, tol)   && lift.all.n   == 12, 'T1 all');
npass = npass + 1;

% touchdown(stance onset) of posterior -> liftoff of anterior (default 'signed').
% stance onset = swing onset + swingdur (=5); nearest RM swing onset within
% +/-0.5P is +5 ms ahead. signed matches all 4 ref onsets -> n=16 for 'all'.
td = computeP2ALag(bd, walk_t0, walk_t1, timestamps, fly, 'stance');
assert(approx(td.RH_to_RM.mean, 5, tol), 'T1 RH_to_RM touchdown lag');
assert(approx(td.all.mean, 5, tol) && td.all.n == 16, 'T1 touchdown all');
npass = npass + 1;

%% Test 2: zero lag (synchronous liftoffs) -> lag 0
onsets = cell(1,6);
for l = 1:6, onsets{l} = 10 + cyc; end
bd = makeboutdata(onsets, swingdur);
lift = computeP2ALag(bd, walk_t0, walk_t1, timestamps, fly, 'swing');
assert(approx(lift.all.mean, 0, tol) && lift.all.n == 12, 'T2 synchronous lag 0');
npass = npass + 1;

%% Test 3: forward-window on a "reversed" wave (anterior leads by 10)
% RM lifts off 10 frames BEFORE RH. Forward-only window [ref,ref+period)
% cannot see the leading onset, so it matches the next cycle: lag = 50-10 = 40.
onsets = cell(1,6);
onsets{3} = 10 + cyc;  onsets{2} = 0 + cyc;  onsets{1} = 30 + cyc;
onsets{4} = 15 + cyc;  onsets{5} = 25 + cyc;  onsets{6} = 35 + cyc;
bd = makeboutdata(onsets, swingdur);
lift = computeP2ALag(bd, walk_t0, walk_t1, timestamps, fly, 'swing');
assert(approx(lift.RH_to_RM.mean, 40, tol), 'T3 reversed -> 40ms (forward window)');
assert(all(lift.RH_to_RM.data >= 0), 'T3 lags non-negative');
npass = npass + 1;

%% Test 4: single reference onset -> NaN mean, n 0
onsets = cell(1,6);
for l = 1:6, onsets{l} = 10 + cyc; end
onsets{3} = 10;   % RH single onset
bd = makeboutdata(onsets, swingdur);
lift = computeP2ALag(bd, walk_t0, walk_t1, timestamps, fly, 'swing');
assert(isnan(lift.RH_to_RM.mean) && lift.RH_to_RM.n == 0, 'T4 single onset -> NaN/0');
npass = npass + 1;

%% Test 5: anterior never swings in window -> NaN for that pair
onsets = cell(1,6);
for l = 1:6, onsets{l} = 10 + cyc; end
onsets{2} = 480;  % RM only swings far outside RH windows
bd = makeboutdata(onsets, swingdur);
lift = computeP2ALag(bd, walk_t0, walk_t1, timestamps, fly, 'swing');
assert(lift.RH_to_RM.n == 0 && isnan(lift.RH_to_RM.mean), 'T5 missing match -> NaN/0');
npass = npass + 1;

%% Test 6: signed pairing returns a NEGATIVE lag when the anterior leg leads
% RH stance onsets [50,100,150] (period 50, half-window 25); RM swing onsets
% [45,105,155]. Nearest to the first ref (50) is 45 -> anterior leads -> -5 ms;
% the next two match +5 ms. Exercises the signed branch + sign convention.
onsets = cell(1,6);
for l = 1:6, onsets{l} = 10 + cyc; end
onsets{3} = [45 95 145];    % RH swing -> stance onsets [50 100 150]
onsets{2} = [45 105 155];   % RM swing onsets
bd = makeboutdata(onsets, swingdur);
td = computeP2ALag(bd, walk_t0, walk_t1, timestamps, fly, 'stance', 'signed');
assert(approx(td.RH_to_RM.data(1), -5, tol), 'T6 signed anterior-leads -> negative');
assert(approx(td.RH_to_RM.data(2),  5, tol) && approx(td.RH_to_RM.data(3), 5, tol), 'T6 signed positive');
npass = npass + 1;

fprintf('test_P2ALag: ALL %d test groups PASSED\n', npass);

end

% ---------------------------------------------------------------------------
function bd = makeboutdata(onsets, swingdur)
% Build currflyboutdata-like struct: .perlimb(l).swing/.stance with
% .start_indices/.end_indices. swing = [onset, onset+swingdur]; stance fills
% the gap to the next onset (touchdown = swing end).
nlimb = numel(onsets);
for l = 1:nlimb
    s = sort(onsets{l}(:)');
    bd.perlimb(l).swing.start_indices = s;
    bd.perlimb(l).swing.end_indices   = s + swingdur;
    % stance onset (touchdown) at swing end; end just before next swing onset
    st0 = s + swingdur;
    if numel(s) >= 2
        st1 = [s(2:end) - 1, st0(end) + 40];
    else
        st1 = st0 + 40;
    end
    bd.perlimb(l).stance.start_indices = st0;
    bd.perlimb(l).stance.end_indices   = st1;
end
end

function tf = approx(a, b, tol)
tf = abs(a - b) <= tol;
end
