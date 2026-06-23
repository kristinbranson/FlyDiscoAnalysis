function perexpgait = computePerExpgaitclass(perwalk_metrics, perfly_metrics)
% computePerExpgaitclass - Aggregate per-walk gait_class structs into per-experiment.
%
% Frame-pooled aggregation: counts and raw class codes summed/concatenated
% across all walks (and equivalently across all flies) in the experiment.
% Speed-conditional pooling is a downstream addition (see plans.md).
%
% Inputs:
%   perwalk_metrics - struct array (one entry per walk across all flies)
%   perfly_metrics  - struct array (one entry per fly), with per-fly counts
%                     already aggregated; used as a cross-check / fly-mean source.
%
% Output:
%   perexpgait struct with:
%     .frm_data_exp        - concatenated class codes across all walks (uint8)
%     .frm_n_exp           - total frame count
%     .tripod_count_exp    - sum across walks
%     .tetrapod_count_exp
%     .grounded_count_exp
%     .airborne_count_exp
%     .other_count_exp
%
% Note: perfly_metrics is currently unused — kept as a parameter to mirror
% computePerExpwalkperframefeatures' signature and to provide a hook for
% per-fly fraction stats if needed later.

class_names = {'tripod', 'tetrapod', 'grounded', 'airborne', 'other'};

nwalks = numel(perwalk_metrics);

% Concatenate raw per-frame class codes across all walks
data_cells = cell(1, nwalks);
for w = 1:nwalks
    data_cells{w} = perwalk_metrics(w).gait_class.data;
end
perexpgait.frm_data_exp = horzcat(data_cells{:});
perexpgait.frm_n_exp    = numel(perexpgait.frm_data_exp);

% Co-indexed smoothed speed (per-frame) for speed-binning; gait codes carry no
% NaN, so frm_speed_exp is the smoothvelmag concatenation in the same frame
% order. [] if unavailable or length mismatch (binning then skipped).
perexpgait.frm_speed_exp = [];
if isfield(perwalk_metrics, 'smoothvelmag_ctr')
    spd_cells = cell(1, nwalks);
    for w = 1:nwalks
        spd_cells{w} = perwalk_metrics(w).smoothvelmag_ctr.data;
    end
    spd = horzcat(spd_cells{:});
    if numel(spd) == numel(perexpgait.frm_data_exp)
        perexpgait.frm_speed_exp = spd;
    end
end

% Sum per-class counts across walks
for c = 1:numel(class_names)
    cf = [class_names{c} '_count'];
    counts = zeros(1, nwalks);
    for w = 1:nwalks
        counts(w) = perwalk_metrics(w).gait_class.(cf);
    end
    perexpgait.([class_names{c} '_count_exp']) = sum(counts);
end

end
