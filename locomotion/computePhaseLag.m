function phaseoffsetdata = computePhaseLag(currflyboutdata, walk_t0, walk_t1)


% input leg x data matrix
phaseoffsetdata = struct;
legorder = {'RF','RM','RH','LH','LM','LF'};
pre_pad_proportion = 0;

% first leg in reference
diffs = [5, 1
    5, 2
    5, 3
    5, 4
    5, 6
    6,1
    5,2
    4,3
    6,5
    1,2
    2,3
    1,5
    6,2
    2,4];

% 5 vs reference limb

for d = 1:numel(diffs)/2
    ref_limb = diffs(d,1);
    % find swing starts in this walk bout
    idx = find(currflyboutdata.perlimb(ref_limb).swing.start_indices>= walk_t0 & currflyboutdata.perlimb(ref_limb).swing.end_indices<= walk_t1);
    refswings_t0s = currflyboutdata.perlimb(ref_limb).swing.start_indices(idx);
    refswings_t1s = currflyboutdata.perlimb(ref_limb).swing.end_indices(idx);

    % check pairing is correct
    assert(all(refswings_t1s - refswings_t0s > 0), 'Invalid swing durations');

    stepdata = cell(1,numel(walk_t0));
    for l = 1:numel(currflyboutdata.perlimb)
        stepdata{l} = currflyboutdata.perlimb(l).swing.start_indices';
    end
    matched_stepdata = findClosestSteps(refswings_t0s', stepdata, pre_pad_proportion);

    periods = diff(refswings_t0s');

    % phaselagdiffs = mod((matched_stepdata(diffs(d,2),:)) - (matched_stepdata(diffs(d,1),:)) ./periods,1.0);
    frame_diff = matched_stepdata(diffs(d,2),:) - matched_stepdata(diffs(d,1),:);

     % Convert to phase in radians using the period
    phase_lag_rad = 2*pi * frame_diff ./ periods;

      % Wrap to [-π, π] using circular distance      
    phaselagdiffs = angle(exp(1i * phase_lag_rad));

    name = [legorder{diffs(d,1)},'_',legorder{diffs(d,2)}];
    currdata = phaselagdiffs;
    phaseoffsetdata.(name).data = currdata;
    phaseoffsetdata.(name).mean = circ_mean(currdata(~isnan(currdata)));
    phaseoffsetdata.(name).std = circ_std(currdata(~isnan(currdata)));
    phaseoffsetdata.(name).n = numel(currdata(~isnan(currdata)));

end
end


