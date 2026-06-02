%% test_LR_symmetrization.m
% Synthetic test of the L/R x-pooling fix in computeboutmetrics2.m.
%
% Constructs synthetic tip_pos_body for 6 limbs across 3 walk bouts with
% distinct patterns:
%   Bout A: medial sweep   - legs converge toward midline during stance
%   Bout B: pure backward  - no sideways component, larger amplitude
%   Bout C: lateral sweep  - legs diverge from midline during stance
%
% Limb convention: limbs 1-3 (RF, RM, RH) sit at +x; limbs 4-6 (LH, LM, LF)
% at -x. L and R within each bout have mirror x trajectories.
%
% Approach:
%   1. Build synthetic tip_pos_body + minimal trx/aptdata
%   2. Run computeStepFeatures per limb -> bout_metrics.perfly(fly).perlimb
%   3. Hand-check a few per-limb step features (AEP_BL, step_direction,
%      amplitude_BL) against expected values
%   4. Replicate the L+R pooling block from computeboutmetrics2.m and build
%      pairs / all_limbs entries
%   5. Assert pair/all entries match the |x|-pooled expectation, per-limb
%      data was not mutated, and non-x fields were not touched by the fix
%
% NOTE: The pooling loop here is a parallel copy of the production code at
% computeboutmetrics2.m lines ~142-180. If that loop is modified, this test
% must be updated in lockstep.

addpath('/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis');
modpath;

fprintf('=== test_LR_symmetrization ===\n\n');

%% ---- 1. Synthetic data ----------------------------------------------------
nlimbs   = 6;
nframes  = 40;
nlandmrk = 6;
fly      = 1;
BL       = 2;   % body length in px: meanbodylength = 4*mean(trx(fly).a) = 4*0.5 = 2

% Stance boundaries (same across all limbs for test simplicity).
% Convention: stance_t1s is EXCLUSIVE (first frame after stance).
stance_t0s = [5; 17; 30];
stance_t1s = [13; 26; 38];     % last stance frames are 12, 25, 37
step_t0s   = [5; 17; 30];
step_t1s   = [17; 30; 38];     % step ends at next AEP (exclusive)
nbouts     = numel(stance_t0s);

% Per-bout R-side AEP/PEP x. L-side is the negation.
%   Bout A: medial sweep (AEPx > PEPx for R, AEPx < PEPx for L)
%   Bout B: pure backward, no sideways
%   Bout C: lateral sweep (AEPx < PEPx for R, AEPx > PEPx for L)
rside_AEP_x = [+2, +2, +1];
rside_PEP_x = [+1, +2, +2];

% Y trajectories (same sign for L and R since A-P axis is shared).
all_AEP_y = [-3, -4, -3];
all_PEP_y = [+3, +4, +3];

% Build tip_pos_body. Linear interpolation between AEP and PEP within stance;
% NaN elsewhere (swing frames are never queried by computeStepFeatures).
tip_pos_body = nan(nlimbs, 2, nframes);
for L = 1:nlimbs
    side_sign = 1 - 2*(L >= 4);          % R-side (L=1-3) -> +1; L-side (L=4-6) -> -1
    for b = 1:nbouts
        t0 = stance_t0s(b);
        t1 = stance_t1s(b) - 1;          % last stance frame
        nstance = t1 - t0 + 1;
        ax = side_sign * rside_AEP_x(b);
        px = side_sign * rside_PEP_x(b);
        ay = all_AEP_y(b);
        py = all_PEP_y(b);
        for k = 0:(nstance-1)
            frac = k / max(nstance - 1, 1);
            tip_pos_body(L, 1, t0 + k) = ax + frac*(px - ax);
            tip_pos_body(L, 2, t0 + k) = ay + frac*(py - ay);
        end
    end
end

% Minimal trx and aptdata so computeStepFeatures runs (see test_PEP_indexing_fix.m).
trx = struct();
trx(fly).a   = ones(1, nframes) * 0.5;   % body half-axis; meanbodylength = 4*mean(a) = 2
trx(fly).off = 0;

aptdata = struct();
aptdata.pTrk = cell(1, fly);
aptdata.pTrk{fly} = nan(nlandmrk, 2, nframes);
for t = 1:nframes
    aptdata.pTrk{fly}(:, 1, t) = t;       % values arbitrary, only length_px uses them
    aptdata.pTrk{fly}(:, 2, t) = t*10;
end
legtip_landmarknums = 1:nlandmrk;
currfly_timestamps  = (0:nframes-1)' / 150;   % 150 fps

%% ---- 2. computeStepFeatures per limb --------------------------------------
bout_metrics = struct();
for L = 1:nlimbs
    sf = computeStepFeatures(fly, trx, aptdata, tip_pos_body, ...
        legtip_landmarknums, L, step_t0s, step_t1s, stance_t0s, stance_t1s, ...
        currfly_timestamps);
    bout_metrics.perfly(fly).perlimb(L).step.stepfeatures = sf;
end

%% ---- 3. Sanity-check per-limb step features -------------------------------

% AEP_BL row 1: signed x. R legs positive, L legs negative.
for L = 1:3
    aex = bout_metrics.perfly(fly).perlimb(L).step.stepfeatures.AEP_BL(1,:);
    expected = rside_AEP_x / BL;
    assert(all(abs(aex - expected) < 1e-12), ...
        'Limb %d (R) AEP_BL(1,:) wrong: got %s, expected %s', ...
        L, mat2str(aex), mat2str(expected));
end
for L = 4:6
    aex = bout_metrics.perfly(fly).perlimb(L).step.stepfeatures.AEP_BL(1,:);
    expected = -rside_AEP_x / BL;
    assert(all(abs(aex - expected) < 1e-12), ...
        'Limb %d (L) AEP_BL(1,:) wrong: got %s, expected %s', ...
        L, mat2str(aex), mat2str(expected));
end
fprintf('Per-limb AEP_BL row 1 signs OK (R+, L-).\n');

% AEP_BL row 2: signed y, same for L and R.
for L = 1:nlimbs
    aey = bout_metrics.perfly(fly).perlimb(L).step.stepfeatures.AEP_BL(2,:);
    expected = all_AEP_y / BL;
    assert(all(abs(aey - expected) < 1e-12), ...
        'Limb %d AEP_BL(2,:) wrong: got %s, expected %s', ...
        L, mat2str(aey), mat2str(expected));
end
fprintf('Per-limb AEP_BL row 2 (signed y) OK.\n');

% step_direction per limb: atan2(Δy, Δx).
% R legs: Δx = rside_PEP_x - rside_AEP_x; L legs: Δx negated -> mirror angle.
for L = 1:3
    sd       = bout_metrics.perfly(fly).perlimb(L).step.stepfeatures.step_direction;
    expected = atan2(all_PEP_y - all_AEP_y, rside_PEP_x - rside_AEP_x);
    assert(all(abs(sd - expected) < 1e-12), ...
        'Limb %d (R) step_direction wrong: got %s, expected %s', ...
        L, mat2str(sd), mat2str(expected));
end
for L = 4:6
    sd       = bout_metrics.perfly(fly).perlimb(L).step.stepfeatures.step_direction;
    expected = atan2(all_PEP_y - all_AEP_y, -(rside_PEP_x - rside_AEP_x));
    assert(all(abs(sd - expected) < 1e-12), ...
        'Limb %d (L) step_direction wrong: got %s, expected %s', ...
        L, mat2str(sd), mat2str(expected));
end
fprintf('Per-limb step_direction OK (R and mirrored L).\n');

% amplitude_BL (Euclidean): sign-invariant; same for L and R.
for L = 1:nlimbs
    amp = bout_metrics.perfly(fly).perlimb(L).step.stepfeatures.amplitude_BL;
    dx  = rside_PEP_x - rside_AEP_x;
    dy  = all_PEP_y - all_AEP_y;
    expected = hypot(dx, dy) / BL;
    assert(all(abs(amp - expected) < 1e-12), ...
        'Limb %d amplitude_BL wrong: got %s, expected %s', ...
        L, mat2str(amp), mat2str(expected));
end
fprintf('Per-limb amplitude_BL OK.\n\n');

%% ---- 4. Replicate the L+R pooling block from computeboutmetrics2 ----------
% Mirrors computeboutmetrics2.m lines ~142-180. Keep in sync.

ml_xy_features = {'AEP', 'AEP_BL', 'PEP', 'PEP_BL'};
pairs          = [1,6; 2,5; 3,4];
state_name     = 'step';

flds = fields(bout_metrics.perfly(fly).perlimb(1).(state_name));
for fld = 1:numel(flds)
    fname = flds{fld};

    % pairs
    for p = 1:size(pairs,1)
        L1 = pairs(p, 1);
        L2 = pairs(p, 2);
        a_obj = bout_metrics.perfly(fly).perlimb(L1).(state_name).(fname);
        b_obj = bout_metrics.perfly(fly).perlimb(L2).(state_name).(fname);
        if isstruct(a_obj)
            subflds = fields(a_obj);
            for sbf = 1:numel(subflds)
                a = a_obj.(subflds{sbf});
                b = b_obj.(subflds{sbf});
                if ismember(subflds{sbf}, ml_xy_features) && size(a,1) == 2
                    a(1,:) = abs(a(1,:));
                    b(1,:) = abs(b(1,:));
                end
                bout_metrics.perfly(fly).pairs(p).(state_name).(fname).(subflds{sbf}) = [a, b];
            end
        else
            bout_metrics.perfly(fly).pairs(p).(state_name).(fname) = [a_obj, b_obj];
        end
    end

    % all_limbs
    a_first = bout_metrics.perfly(fly).perlimb(1).(state_name).(fname);
    if isstruct(a_first)
        subflds = fields(a_first);
        for sbf = 1:numel(subflds)
            is_ml = ismember(subflds{sbf}, ml_xy_features);
            currfly_all = cell(1, nlimbs);
            for L = 1:nlimbs
                curr = bout_metrics.perfly(fly).perlimb(L).(state_name).(fname).(subflds{sbf});
                if is_ml && size(curr,1) == 2
                    curr(1,:) = abs(curr(1,:));
                end
                currfly_all{L} = curr;
            end
            bout_metrics.perfly(fly).all_limbs.(state_name).(fname).(subflds{sbf}) = horzcat(currfly_all{:});
        end
    else
        currfly_all = cell(1, nlimbs);
        for L = 1:nlimbs
            currfly_all{L} = bout_metrics.perfly(fly).perlimb(L).(state_name).(fname);
        end
        bout_metrics.perfly(fly).all_limbs.(state_name).(fname) = horzcat(currfly_all{:});
    end
end

%% ---- 5. Pair-level assertions ---------------------------------------------

% AEP_BL row 1: |x| from R then L (both equal |x| since they're mirrors).
expected_AEP_BLx_pair = [rside_AEP_x / BL, rside_AEP_x / BL];
for p = 1:size(pairs,1)
    row1 = bout_metrics.perfly(fly).pairs(p).step.stepfeatures.AEP_BL(1,:);
    assert(all(abs(row1 - expected_AEP_BLx_pair) < 1e-12), ...
        'Pair %d AEP_BL(1,:) (|x|) wrong: got %s, expected %s', ...
        p, mat2str(row1), mat2str(expected_AEP_BLx_pair));
end
fprintf('Pair AEP_BL row 1 (|x|) OK.\n');

% AEP_BL row 2: signed y, untouched.
expected_AEP_BLy_pair = [all_AEP_y / BL, all_AEP_y / BL];
for p = 1:size(pairs,1)
    row2 = bout_metrics.perfly(fly).pairs(p).step.stepfeatures.AEP_BL(2,:);
    assert(all(abs(row2 - expected_AEP_BLy_pair) < 1e-12), ...
        'Pair %d AEP_BL(2,:) (signed y) wrong: got %s, expected %s', ...
        p, mat2str(row2), mat2str(expected_AEP_BLy_pair));
end
fprintf('Pair AEP_BL row 2 (signed y) OK.\n');

% PEP_BL: same pattern.
expected_PEP_BLx_pair = [rside_PEP_x / BL, rside_PEP_x / BL];
expected_PEP_BLy_pair = [all_PEP_y / BL, all_PEP_y / BL];
for p = 1:size(pairs,1)
    row1 = bout_metrics.perfly(fly).pairs(p).step.stepfeatures.PEP_BL(1,:);
    row2 = bout_metrics.perfly(fly).pairs(p).step.stepfeatures.PEP_BL(2,:);
    assert(all(abs(row1 - expected_PEP_BLx_pair) < 1e-12), 'Pair %d PEP_BL(1,:) wrong', p);
    assert(all(abs(row2 - expected_PEP_BLy_pair) < 1e-12), 'Pair %d PEP_BL(2,:) wrong', p);
end
fprintf('Pair PEP_BL rows 1 (|x|) and 2 (y) OK.\n');

% step_direction (1-row): not in ml_xy_features, so plain concat of signed
% per-limb values; in particular L and R values are NOT made equal by the fix.
for p = 1:size(pairs,1)
    L1 = pairs(p, 1);
    L2 = pairs(p, 2);
    expected_sd = [bout_metrics.perfly(fly).perlimb(L1).step.stepfeatures.step_direction, ...
                   bout_metrics.perfly(fly).perlimb(L2).step.stepfeatures.step_direction];
    got_sd      = bout_metrics.perfly(fly).pairs(p).step.stepfeatures.step_direction;
    assert(all(abs(got_sd - expected_sd) < 1e-12), ...
        'Pair %d step_direction wrong (should be signed concat, unchanged by fix)', p);
end
fprintf('Pair step_direction (signed, unchanged by fix) OK.\n');

% amplitude_BL (1-row, sign-invariant): unchanged.
for p = 1:size(pairs,1)
    L1 = pairs(p, 1);
    L2 = pairs(p, 2);
    expected_amp = [bout_metrics.perfly(fly).perlimb(L1).step.stepfeatures.amplitude_BL, ...
                    bout_metrics.perfly(fly).perlimb(L2).step.stepfeatures.amplitude_BL];
    got_amp      = bout_metrics.perfly(fly).pairs(p).step.stepfeatures.amplitude_BL;
    assert(all(abs(got_amp - expected_amp) < 1e-12), ...
        'Pair %d amplitude_BL wrong', p);
end
fprintf('Pair amplitude_BL (unchanged) OK.\n\n');

%% ---- 6. all_limbs assertions ----------------------------------------------

% AEP_BL row 1: |x| concat of all 6 limbs; same |x| in our design.
expected_AEP_BLx_all = repmat(rside_AEP_x / BL, 1, nlimbs);
all_row1 = bout_metrics.perfly(fly).all_limbs.step.stepfeatures.AEP_BL(1,:);
assert(all(abs(all_row1 - expected_AEP_BLx_all) < 1e-12), ...
    'all_limbs AEP_BL(1,:) wrong: got %s, expected %s', ...
    mat2str(all_row1), mat2str(expected_AEP_BLx_all));

% Row 2: signed y, untouched.
expected_AEP_BLy_all = repmat(all_AEP_y / BL, 1, nlimbs);
all_row2 = bout_metrics.perfly(fly).all_limbs.step.stepfeatures.AEP_BL(2,:);
assert(all(abs(all_row2 - expected_AEP_BLy_all) < 1e-12), ...
    'all_limbs AEP_BL(2,:) wrong: got %s, expected %s', ...
    mat2str(all_row2), mat2str(expected_AEP_BLy_all));
fprintf('all_limbs AEP_BL rows 1 (|x|) and 2 (signed y) OK.\n\n');

%% ---- 7. Per-limb data must be unchanged by the pooling --------------------
for L = 1:3
    row1 = bout_metrics.perfly(fly).perlimb(L).step.stepfeatures.AEP_BL(1,:);
    assert(all(row1 > 0), 'Per-limb R AEP_BL row 1 mutated! got %s', mat2str(row1));
end
for L = 4:6
    row1 = bout_metrics.perfly(fly).perlimb(L).step.stepfeatures.AEP_BL(1,:);
    assert(all(row1 < 0), 'Per-limb L AEP_BL row 1 mutated! got %s', mat2str(row1));
end
fprintf('Per-limb data not mutated.\n\n');

fprintf('=== ALL TESTS PASSED ===\n');
