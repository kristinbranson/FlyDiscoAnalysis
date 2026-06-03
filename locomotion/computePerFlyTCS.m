function perflytcs = computePerFlyTCS(walkfeaturestruct)
% computePerFlyTCS - Aggregate per-walk TCS structs into per-fly.
%
% Input:
%   walkfeaturestruct - struct array, each element has .TCS with .both
%                       (combined LM+RM tripod aggregation), which has
%                       .data (1 x n event TCS values in [0,1]),
%                       .event_speeds, .n_steps, .n_nontripod.
%
% Output:
%   perflytcs struct with:
%     .data          - concatenated TCS event values across all walks
%     .event_speeds  - concatenated per-event speeds across all walks
%     .mean, .std, .n - aggregate over events (linear, TCS is scalar [0,1])
%     .n_steps       - summed candidate steps across walks
%     .n_nontripod   - summed non-tripod candidate steps across walks
%
% TCS is a scalar in [0,1] so linear (not circular) statistics are used.

nwalks = numel(walkfeaturestruct);

data_cells   = cell(1, nwalks);
speed_cells  = cell(1, nwalks);
n_steps      = 0;
n_nontripod  = 0;
for w = 1:nwalks
    both = walkfeaturestruct(w).TCS.both;
    data_cells{w}  = both.data;
    speed_cells{w} = both.event_speeds;
    n_steps        = n_steps + both.n_steps;
    n_nontripod    = n_nontripod + both.n_nontripod;
end

perflytcs.data         = horzcat(data_cells{:});
perflytcs.event_speeds = horzcat(speed_cells{:});
perflytcs.mean         = mean(perflytcs.data, 'omitnan');
perflytcs.std          = std(perflytcs.data, 'omitnan');
perflytcs.n            = numel(perflytcs.data);
perflytcs.n_steps      = n_steps;
perflytcs.n_nontripod  = n_nontripod;

end
