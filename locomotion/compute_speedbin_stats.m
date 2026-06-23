function binned = compute_speedbin_stats(values, speeds, edges, stattype, labels)
% COMPUTE_SPEEDBIN_STATS  Per-speed-bin sufficient statistics.
%   Splits co-indexed (values, speeds) into speed bins and returns the same
%   (mean,std,Z) / (frac,Z) summary used elsewhere in computeStatsPerExp, one
%   per bin. Used by both binning paths (per-event and per-frame).
%
%   values   1xN metric values (per event or per frame)
%   speeds   1xN co-indexed speed (mm/s); frames/events with NaN speed are
%            unassigned (dropped from every bin)
%   edges    interior bin edges, e.g. [8 20] -> 3 bins:
%            [-inf,8) [8,20) [20,inf)   (right-open; matches the plots)
%   stattype 'linear'   -> .mean (omitnan), .std (omitnan), .Z (n non-NaN)
%            'circular' -> .mean (circ_mean), .std (circ_std), .Z
%            'fraction' -> values are 0/1 membership; .frac (mean), .Z
%   labels   optional 1x(numel(edges)+1) cellstr of bin names; default
%            {'slow','med','fast'} for 2 edges, else {'bin1',...}.
%
%   Returns struct `binned` with one field per label, each a struct shaped
%   like the unbinned field for that stattype. Empty bins -> NaN stats, Z=0.

if nargin < 5 || isempty(labels)
    nbin = numel(edges) + 1;
    if nbin == 3
        labels = {'slow','med','fast'};
    else
        labels = arrayfun(@(k) sprintf('bin%d',k), 1:nbin, 'uni',0);
    end
end
nbin = numel(labels);
assert(nbin == numel(edges)+1, 'labels must have numel(edges)+1 entries');

values = values(:)';
speeds = speeds(:)';
binidx = discretize(speeds, [-inf, edges, inf]);   % 1..nbin, NaN if speed NaN

binned = struct();
for b = 1:nbin
    sel = (binidx == b);
    v = values(sel);
    v = v(~isnan(v));
    s = struct();
    switch stattype
        case 'linear'
            if isempty(v)
                s.mean = NaN; s.std = NaN; s.Z = 0;
            else
                s.mean = mean(v); s.std = std(v); s.Z = numel(v);
            end
        case 'circular'
            if isempty(v)
                s.mean = NaN; s.std = NaN; s.Z = 0;
            else
                s.mean = circ_mean(v(:)); s.std = circ_std(v(:)); s.Z = numel(v);
            end
        case 'fraction'
            if isempty(v)
                s.frac = NaN; s.Z = 0;
            else
                s.frac = mean(v); s.Z = numel(v);
            end
        otherwise
            error('Unknown stattype "%s". Use linear|circular|fraction.', stattype);
    end
    binned.(labels{b}) = s;
end
