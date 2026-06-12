function out = computeP2ALag(currflyboutdata, walk_t0, walk_t1, timestamps, fly, ref_event_type, pairing)
% computeP2ALag - Per-walk posterior-to-anterior (P2A) onset lag (ms).
%
% Time delay from a posterior (reference) leg event to a nearby swing onset
% (liftoff) of the anterior (target) leg, for the four ipsilateral leg pairs.
%
% ref_event_type selects the reference event on the posterior leg:
%   'swing'  -> posterior liftoff   => Pliftoff2Aliftoff_lag
%   'stance' -> posterior touchdown => Ptouchdown2Aliftoff_lag
% The target event is always the anterior leg's liftoff (swing onset).
%
% pairing selects how the anterior onset is matched:
%   'forward' : nearest onset in the forward window [ref, ref+period); lag >= 0.
%               Right for the liftoff->liftoff wave, which is ~antiphase
%               (~0.5*period): a forward-only window is unimodal there, whereas
%               a signed/nearest rule would split the antiphase population into
%               +0.5 / -0.5 modes (wrap artifact).
%   'signed'  : nearest onset in the symmetric window [ref-0.5P, ref+0.5P],
%               P = median reference period; signed lag (negative = anterior
%               leads). Right for touchdown->liftoff (Cruse Rule 2), which is
%               centered ~0; a forward-only window would push slight-negative
%               events to ~+1 period (spurious large lags).
% Default: 'forward' for ref 'swing', 'signed' for ref 'stance'.
%
% Window choice rationale (empirical, all-floor-walks histograms): put the
% wrap point away from where the data sits. The +/-0.5P cap also makes the
% match effectively 1:1 and monotonic (adjacent cap windows tile without
% overlap), so no global assignment is needed for regular stepping.
%
% Output struct: one field per pair (RH_to_RM, RM_to_RF, LH_to_LM, LM_to_LF)
% plus segment aggregates H_to_M, M_to_F and pooled all. Each has:
%   .data (per-event lags, ms; NaN where unmatched), .mean, .std, .n (omitnan).
% Aggregates concatenate their constituent pairs' .data then take the stats.
%
% See also: computePhaseLag, findClosestSteps, computePerFlyP2ALag

if nargin < 6 || isempty(ref_event_type), ref_event_type = 'swing'; end
if nargin < 7 || isempty(pairing)
    if strcmp(ref_event_type, 'stance'), pairing = 'signed'; else, pairing = 'forward'; end
end
assert(ismember(ref_event_type, {'swing','stance'}), 'ref_event_type must be ''swing'' or ''stance''');
assert(ismember(pairing, {'forward','signed'}), 'pairing must be ''forward'' or ''signed''');

legorder = {'RF','RM','RH','LH','LM','LF'};

% posterior (reference) -> anterior (target) limb indices in legorder
pairs = [3 2    % RH_to_RM
         2 1    % RM_to_RF
         4 5    % LH_to_LM
         5 6];  % LM_to_LF

% target onsets are always swing onsets (liftoff)
nlimbs = numel(currflyboutdata.perlimb);
stepdata = cell(1, nlimbs);
for l = 1:nlimbs
    stepdata{l} = currflyboutdata.perlimb(l).swing.start_indices';
end

nts = numel(timestamps);
out = struct;

for d = 1:size(pairs,1)
    ref_limb    = pairs(d,1);
    target_limb = pairs(d,2);
    name = [legorder{ref_limb}, '_to_', legorder{target_limb}];

    refbouts = currflyboutdata.perlimb(ref_limb).(ref_event_type);
    idx = find(refbouts.start_indices >= walk_t0 & refbouts.end_indices <= walk_t1);
    ref_onsets = sort(refbouts.start_indices(idx));

    if strcmp(pairing, 'forward')
        lag = pair_forward(ref_onsets, stepdata, target_limb, timestamps, nts, ...
                           fly, walk_t0, walk_t1, name, ref_event_type);
    else
        lag = pair_signed(ref_onsets, stepdata{target_limb}, timestamps, nts);
    end
    out.(name) = lagstats(lag);
end

% segment aggregates and pooled all (concatenate-then-mean)
out.H_to_M = lagstats([out.RH_to_RM.data, out.LH_to_LM.data]);
out.M_to_F = lagstats([out.RM_to_RF.data, out.LM_to_LF.data]);
out.all    = lagstats([out.RH_to_RM.data, out.RM_to_RF.data, ...
                       out.LH_to_LM.data, out.LM_to_LF.data]);

end

% ---------------------------------------------------------------------------
function lag = pair_forward(ref_onsets, stepdata, target_limb, ts, nts, fly, walk_t0, walk_t1, name, ref_event_type)
% forward-only [ref, ref+period) match via findClosestSteps; lag >= 0
try
    matched_stepdata = findClosestSteps(ref_onsets', stepdata, 0);   % pre_pad=0
    matched_target = matched_stepdata(target_limb, :);
catch ME
    warning('FlyDisco:P2ALag:pairingError', ...
        'computeP2ALag(%s,forward): findClosestSteps failed for fly %d, frames %d-%d, pair %s: %s', ...
        ref_event_type, fly, walk_t0, walk_t1, name, ME.message);
    matched_target = nan(1, max(0, numel(ref_onsets)-1));
end
matched_ref = ref_onsets(1:end-1); matched_ref = matched_ref(:)';
lag = nan(1, numel(matched_target));
valid = ~isnan(matched_target) & ~isnan(matched_ref) & ...
        matched_target >= 1 & matched_target <= nts & ...
        matched_ref    >= 1 & matched_ref    <= nts;
lag(valid) = (ts(matched_target(valid)) - ts(matched_ref(valid))) * 1000;
end

% ---------------------------------------------------------------------------
function lag = pair_signed(ref_onsets, ant_onsets, ts, nts)
% nearest anterior onset within +/-0.5*median(period); signed lag (ms)
n = numel(ref_onsets);
lag = nan(1, n);
if n < 2, return; end                       % need >=2 onsets for a period
ant = sort(ant_onsets(:)');
halfwin = 0.5 * median(diff(ref_onsets));
for s = 1:n
    r = ref_onsets(s);
    cand = ant(ant >= r - halfwin & ant <= r + halfwin);
    if isempty(cand), continue; end
    [~, k] = min(abs(cand - r));            % nearest (prev or next)
    m = cand(k);
    if m >= 1 && m <= nts && r >= 1 && r <= nts
        lag(s) = (ts(m) - ts(r)) * 1000;    % signed
    end
end
end

% ---------------------------------------------------------------------------
function s = lagstats(data)
% linear stats over per-event lags (ms), ignoring NaN
s.data = data;
s.mean = mean(data, 'omitnan');
s.std  = std(data, 'omitnan');
s.n    = sum(~isnan(data));
end
