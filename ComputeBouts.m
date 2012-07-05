function [bout_starts,bout_ends] = ComputeBouts(isbehavior,doanalyze)

[starts_behavior,ends_behavior] = get_interval_ends(isbehavior);
tsanalyze = find(doanalyze);

% assign each bout to its center frame
mids_behavior = ceil((starts_behavior + ends_behavior-1)/2);

% only care about bouts with centers in tsanalyze
idx = ismember(mids_behavior,tsanalyze);
bout_starts = starts_behavior(idx);
bout_ends = ends_behavior(idx)-1;
