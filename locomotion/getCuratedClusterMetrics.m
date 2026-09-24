function keys = getCuratedClusterMetrics()
% GETCURATEDCLUSTERMETRICS  Curated metric set for line clustering and the
% line-comparison plots (ON-vs-delta, etc). This is getCuratedSpeedbinMetrics()
% MINUS gait-class fractions and nfeet_ground (dropped 2026-06-26: gait fractions
% are compositional/collinear and nfeet is largely redundant with them).
%
% Returns LED-stripped keys 'feature__state__qualifier' (no LED, no bin). The
% clustering/plotting feature space is these keys x {unbinned,slow,med,fast} x
% {LEDon, delta}. 64 base features.

keys = getCuratedSpeedbinMetrics();
drop = startsWith(keys, 'gait_class__') | startsWith(keys, 'nfeet_ground__');
keys = keys(~drop);
assert(numel(keys) == 64, 'expected 64 cluster keys, got %d', numel(keys));
end
