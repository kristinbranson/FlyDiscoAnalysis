function stats = CombinePerExpStats(statsperexp)

% which exps we will look at
goodidx = ~isnan(statsperexp.Z) & ~isnan(statsperexp.meanmean) & ...
  ~isinf(statsperexp.meanmean);

stats = struct;
stats.Z = sum(statsperexp.Z(goodidx));
stats.meanZ = mean(statsperexp.Z(goodidx));
stats.meanmean = mean(statsperexp.meanmean(goodidx));
stats.stdmean = std(statsperexp.meanmean(goodidx),1);
stats.meanstd = mean(statsperexp.stdmean(goodidx));
stats.meanprctiles = mean(statsperexp.meanprctiles(goodidx,:),1);
stats.stdprctiles = std(statsperexp.meanprctiles(goodidx,:),1,1);
stats.n = nnz(goodidx);