function statsperexp = CombinePerFrameStats(statsperfly)

nprctiles = size(statsperfly.prctiles,1);

% which flies we will look at
goodidx = ~isnan(statsperfly.Z) & ~isnan(statsperfly.mean) & ...
  ~isinf(statsperfly.mean);

% total number of frames of data analyzed in this experiment
statsperexp.Z = sum(statsperfly.Z(goodidx));

% average number of frames of data analyzed per fly, weighted by weight of
% fly
statsperexp.meanZ = weighted_mean(statsperfly.Z(goodidx)',...
  statsperfly.fracframesanalyzed(goodidx)');

% mean of per-fly means
[statsperexp.meanmean,s] = weighted_mean_cov(statsperfly.mean(goodidx)',...
  statsperfly.fracframesanalyzed(goodidx)');

% standard deviation of per-fly means
statsperexp.stdmean = sqrt(s);

% mean of per-fly standard deviations
statsperexp.meanstd = weighted_mean(statsperfly.std(goodidx)',statsperfly.fracframesanalyzed(goodidx)');

% mean and standard deviations of per-fly percentiles
statsperexp.meanprctiles = nan(1,nprctiles);
statsperexp.stdprctiles = nan(1,nprctiles);
for j = 1:nprctiles,
  [statsperexp.meanprctiles(j),s] = weighted_mean_cov(statsperfly.prctiles(j,goodidx)',statsperfly.fracframesanalyzed(goodidx)');
  statsperexp.stdprctiles(j) = sqrt(s);
end