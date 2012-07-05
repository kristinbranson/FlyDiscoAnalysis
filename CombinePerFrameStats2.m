function statsperexp = CombinePerFrameStats2(statsperfly,minZboth,minZfly)

nprctiles = size(statsperfly.prctiles,1);

% which flies we will look at
goodidx = ~isnan(statsperfly.Z) & ~isnan(statsperfly.mean) & ...
  ~isinf(statsperfly.mean) & (statsperfly.Z >= minZboth) & (statsperfly.Zfly >= minZfly);

% total number of frames of data analyzed in this experiment
statsperexp.Z = sum(statsperfly.Z(goodidx));

% average number of frames of data analyzed per fly, weighted by weight of
% fly
statsperexp.meanZ = weighted_mean(statsperfly.Z(goodidx)',...
  statsperfly.fracframesanalyzed(goodidx)');

% mean of per-fly means, the total weight for a trajectory is the number of
% frames for which doanalyze_fly is true
statsperexp.meanmean = weighted_mean(statsperfly.mean(goodidx)',...
  statsperfly.Zfly(goodidx)');

% scatter around the meanmean, but the weighting is slightly different --
% the weight for each trajectory is the *fraction* of frames for which
% doanalyze_fly is true
[s] = weighted_scatter(statsperfly.mean(goodidx)',...
  statsperexp.meanmean,statsperfly.fracframesanalyzed(goodidx)');

% standard deviation of per-fly means
statsperexp.stdmean = sqrt(s);

% mean of per-fly standard deviations
statsperexp.meanstd = weighted_mean(statsperfly.std(goodidx)',statsperfly.Zfly(goodidx)');

% mean and standard deviations of per-fly percentiles
statsperexp.meanprctiles = nan(1,nprctiles);
statsperexp.stdprctiles = nan(1,nprctiles);
for j = 1:nprctiles,
  statsperexp.meanprctiles(j) = weighted_mean(statsperfly.prctiles(j,goodidx)',statsperfly.Zfly(goodidx)');
  s = weighted_scatter(statsperfly.prctiles(j,goodidx)',statsperexp.meanprctiles(j),statsperfly.fracframesanalyzed(goodidx)');
  statsperexp.stdprctiles(j) = sqrt(s);
end

statsperexp.n = numel(goodidx);