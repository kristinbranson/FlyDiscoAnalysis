function perexptcs = computePerExpTCS(perwalk_metrics, perfly_metrics)
% computePerExpTCS - Aggregate per-walk TCS structs into per-experiment.
%
% Event-pooled aggregation: TCS event values and per-event speeds are
% concatenated across all walks (equivalently across all flies) in the
% experiment, mirroring computePerExpgaitclass' frame-pooled approach.
%
% Inputs:
%   perwalk_metrics - struct array (one entry per walk across all flies),
%                     each with .TCS.both.
%   perfly_metrics  - struct array (one entry per fly); currently unused,
%                     kept to mirror computePerExpgaitclass' signature and
%                     provide a hook for per-fly TCS stats if needed later.
%
% Output:
%   perexptcs struct with:
%     .data          - concatenated TCS event values across all walks
%     .event_speeds  - concatenated per-event speeds across all walks
%     .mean, .std, .n - aggregate over events (linear, TCS is scalar [0,1])
%     .n_steps       - summed candidate steps across walks
%     .n_nontripod   - summed non-tripod candidate steps across walks

nwalks = numel(perwalk_metrics);

data_cells  = cell(1, nwalks);
speed_cells = cell(1, nwalks);
n_steps     = 0;
n_nontripod = 0;
for w = 1:nwalks
    both = perwalk_metrics(w).TCS.both;
    data_cells{w}  = both.data;
    speed_cells{w} = both.event_speeds;
    n_steps        = n_steps + both.n_steps;
    n_nontripod    = n_nontripod + both.n_nontripod;
end

perexptcs.data         = horzcat(data_cells{:});
perexptcs.event_speeds = horzcat(speed_cells{:});
perexptcs.mean         = mean(perexptcs.data, 'omitnan');
perexptcs.std          = std(perexptcs.data, 'omitnan');
perexptcs.n            = numel(perexptcs.data);
perexptcs.n_steps      = n_steps;
perexptcs.n_nontripod  = n_nontripod;

end
