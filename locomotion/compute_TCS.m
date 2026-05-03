function TCS = compute_TCS(currflyboutdata, walk_t0, walk_t1, loctall, locball, speed_data)
% compute_TCS - Tripod Coordination Strength per Wosnitza et al. 2012
% (J Exp Biol 216(3):480-491). Computes TCS = t2/t1 for tripod swing
% events within a single walking bout.
%
% Definition (Wosnitza p. 483):
%   For each candidate tripod step:
%     t1 = max(swing offsets across the 3 legs) - min(swing onsets across the 3 legs)
%     t2 = max(0, min(offsets) - max(onsets))    % all-3-in-swing duration
%     TCS = t2 / t1                              % in [0, 1]
%   A step qualifies as a tripod event iff t2 >= 1 frame.
%
% Tripod groupings (Strauss 1990):
%   Tripod LM = {LM, RF, RH}  (lead = LM, leg index 5)
%   Tripod RM = {RM, LF, LH}  (lead = RM, leg index 2)
%
% Lead-leg choice: middle legs are used as the reference leg, deviating
% from Wosnitza's "front leg with most cycles" convention. Reason:
% APT front-leg tracking has elevated noise from head overlap; middle legs
% give cleaner step cycles in this setup.
%
% Per-event speed: mean of raw `velmag_ctr` over the t1 envelope
% [min(swing_starts), max(swing_ends)-1]. Wosnitza smoothed with a
% 5-frame gliding average at 500 fps (= 10 ms); the per-event averaging
% over t1 (typically 20-40 ms = 3-6 frames at 150 fps) provides
% comparable smoothing without an upstream filter.
%
% Restriction C (Hilbert-valid window): candidate steps must have their
% start_index within the intersection of all 3 tripod legs'
% [first peak, last peak] windows, mirroring the restriction
% computeContinuousPhaseDiff_hilbert.m applies to phase data.
%
% Inputs:
%   currflyboutdata - obj.limbBoutData(fly); .perlimb(L).{swing,stance,step}.{start,end}_indices
%                     Indices in per-fly trajectory format. End indices are
%                     EXCLUSIVE (first frame after the bout) per detect_bouts.
%   walk_t0, walk_t1 - trajectory-frame indices for this walking bout
%   loctall  - {1 x nlimb} cell, each [peak_value; local_frame_idx] of top peaks
%              of normalized Y tip position (local indices = within walk window)
%   locball  - {1 x nlimb} cell, same for bottom peaks
%   speed_data - per-trajectory-frame speed values (e.g. velmag_ctr), full length
%
% Output: TCS struct with fields tripod_LM, tripod_RM, both. Each contains:
%   .data         - per-event TCS values (1 x n)
%   .event_times  - trajectory frame of step start (= lead-leg AEP touchdown)
%   .event_speeds - mean speed over t1 envelope per event
%   .step_idx     - index into perlimb(lead).step.start_indices for each event
%   .mean, .std, .n - aggregate stats over the kept events
%   .n_steps      - total candidate steps considered (in valid window)
%   .n_nontripod  - candidate steps that did NOT qualify as tripod events
%
% .both also includes:
%   .tripod_id    - 1 for LM-led events, 2 for RM-led events
%
% Note: tripod fraction = n / n_steps (= 1 - n_nontripod/n_steps).

% leg indices in tip_pos_body convention: RF=1, RM=2, RH=3, LH=4, LM=5, LF=6
LM_tripod = struct('lead', 5, 'others', [1, 3]);   % LM, RF, RH
RM_tripod = struct('lead', 2, 'others', [6, 4]);   % RM, LF, LH

TCS.tripod_LM = compute_TCS_one_group(currflyboutdata, walk_t0, walk_t1, ...
                                      loctall, locball, speed_data, LM_tripod);
TCS.tripod_RM = compute_TCS_one_group(currflyboutdata, walk_t0, walk_t1, ...
                                      loctall, locball, speed_data, RM_tripod);

% Combined aggregation across both tripods
TCS.both.data         = [TCS.tripod_LM.data,         TCS.tripod_RM.data];
TCS.both.event_times  = [TCS.tripod_LM.event_times,  TCS.tripod_RM.event_times];
TCS.both.event_speeds = [TCS.tripod_LM.event_speeds, TCS.tripod_RM.event_speeds];
TCS.both.step_idx     = [TCS.tripod_LM.step_idx,     TCS.tripod_RM.step_idx];
TCS.both.tripod_id    = [ones(1, numel(TCS.tripod_LM.data)), ...
                         2*ones(1, numel(TCS.tripod_RM.data))];
TCS.both.mean         = mean(TCS.both.data, 'omitnan');
TCS.both.std          = std(TCS.both.data, 'omitnan');
TCS.both.n            = numel(TCS.both.data);
TCS.both.n_steps      = TCS.tripod_LM.n_steps + TCS.tripod_RM.n_steps;
TCS.both.n_nontripod  = TCS.tripod_LM.n_nontripod + TCS.tripod_RM.n_nontripod;

end


function tcs_one = compute_TCS_one_group(currflyboutdata, walk_t0, walk_t1, ...
                                         loctall, locball, speed_data, tripod_def)
% Compute TCS for one tripod group, anchored on its middle/lead leg's step cycles.

lead = tripod_def.lead;
others = tripod_def.others;
all_three = [lead, others];

% --- Restriction C: per-leg Hilbert-valid window in trajectory frame ---
sidx_traj = nan(1, 3);
eidx_traj = nan(1, 3);
for j = 1:3
    L = all_three(j);
    if isempty(loctall{L}) || isempty(locball{L})
        tcs_one = empty_tcs_struct();
        return;
    end
    sidx_local = min(loctall{L}(2,1), locball{L}(2,1));
    eidx_local = max(loctall{L}(2,end), locball{L}(2,end));
    sidx_traj(j) = sidx_local + walk_t0 - 1;
    eidx_traj(j) = eidx_local + walk_t0 - 1;
end
valid_window_start = max(sidx_traj);
valid_window_end   = min(eidx_traj);

% --- Lead leg's step cycles ---
step_starts  = currflyboutdata.perlimb(lead).step.start_indices(:);
step_ends    = currflyboutdata.perlimb(lead).step.end_indices(:);
swing_starts_lead = currflyboutdata.perlimb(lead).swing.start_indices(:);
swing_ends_lead   = currflyboutdata.perlimb(lead).swing.end_indices(:);

% Filter steps:
%   - step start within the walk window (matches existing convention)
%   - lead-leg swing onset within Restriction C window (the swing is the
%     meaningful event time for TCS; touchdown can precede the first
%     Hilbert peak without making the swing region unreliable)
in_window = (step_starts >= walk_t0) & (step_starts < walk_t1) & ...
            (swing_starts_lead >= valid_window_start) & ...
            (swing_starts_lead <= valid_window_end);
candidate_idx = find(in_window);
n_candidates = numel(candidate_idx);

% --- Iterate over candidate steps ---
tcs_data     = nan(1, n_candidates);
event_times  = nan(1, n_candidates);
event_speeds = nan(1, n_candidates);
step_idx_out = nan(1, n_candidates);
n_nontripod  = 0;

% Pre-fetch other-leg swing arrays (avoids repeated struct dereference)
other_swing = cell(1, 2);
for j = 1:2
    L = others(j);
    other_swing{j}.starts = currflyboutdata.perlimb(L).swing.start_indices(:);
    other_swing{j}.ends   = currflyboutdata.perlimb(L).swing.end_indices(:);
end

for c = 1:n_candidates
    i = candidate_idx(c);
    step_start = step_starts(i);
    step_end   = step_ends(i);

    swing_starts = nan(3, 1);
    swing_ends   = nan(3, 1);

    % Lead leg's swing within this step (same indexing as step in
    % limbSwingStanceStep: swing[i] is the swing inside step[i])
    swing_starts(1) = currflyboutdata.perlimb(lead).swing.start_indices(i);
    swing_ends(1)   = currflyboutdata.perlimb(lead).swing.end_indices(i);

    valid = true;
    for j = 1:2
        starts_k = other_swing{j}.starts;
        ends_k   = other_swing{j}.ends;

        % Find swings whose start is within step window
        in_step = (starts_k >= step_start) & (starts_k < step_end);
        in_step_idx = find(in_step);
        if isempty(in_step_idx)
            valid = false;
            break;
        elseif numel(in_step_idx) == 1
            idx = in_step_idx;
        else
            % Multiple swings in window: pick max overlap with lead's swing
            cs = starts_k(in_step_idx);
            ce = ends_k(in_step_idx);
            overlap = max(0, min(ce, swing_ends(1)) - max(cs, swing_starts(1)));
            [~, maxi] = max(overlap);
            idx = in_step_idx(maxi);
        end
        swing_starts(j+1) = starts_k(idx);
        swing_ends(j+1)   = ends_k(idx);
    end

    if ~valid
        n_nontripod = n_nontripod + 1;
        continue;
    end

    t1 = max(swing_ends) - min(swing_starts);
    t2 = max(0, min(swing_ends) - max(swing_starts));

    if t2 < 1
        n_nontripod = n_nontripod + 1;
        continue;
    end

    % Valid tripod event
    tcs_data(c)     = t2 / t1;
    event_times(c)  = step_start;
    step_idx_out(c) = i;

    % Per-event speed = mean velmag_ctr over t1 envelope.
    % swing_ends are exclusive, so envelope is [min_start, max_end - 1] inclusive.
    env_start = min(swing_starts);
    env_end   = max(swing_ends) - 1;
    if env_start <= env_end && env_end <= numel(speed_data)
        event_speeds(c) = mean(speed_data(env_start:env_end), 'omitnan');
    end
end

% Keep only valid events in the per-event arrays
keep = ~isnan(tcs_data);
tcs_one.data         = tcs_data(keep);
tcs_one.event_times  = event_times(keep);
tcs_one.event_speeds = event_speeds(keep);
tcs_one.step_idx     = step_idx_out(keep);
tcs_one.mean         = mean(tcs_one.data, 'omitnan');
tcs_one.std          = std(tcs_one.data, 'omitnan');
tcs_one.n            = numel(tcs_one.data);
tcs_one.n_steps      = n_candidates;
tcs_one.n_nontripod  = n_nontripod;

end


function s = empty_tcs_struct()
s.data         = [];
s.event_times  = [];
s.event_speeds = [];
s.step_idx     = [];
s.mean         = NaN;
s.std          = NaN;
s.n            = 0;
s.n_steps      = 0;
s.n_nontripod  = 0;
end
