function histperexp = CombinePerFrameHists(histperfly)

% which flies we will look at
goodidx = ~isnan(histperfly.Z);
% total number of frames of data analyzed in this experiment
histperexp.Z = sum(histperfly.Z(goodidx));
% average number of frames of data analyzed per fly, weighted by weight of
% fly
histperexp.meanZ = weighted_mean(histperfly.Z(goodidx)',...
  histperfly.fracframesanalyzed(goodidx)');
% mean, std fractions
nbins_linear = size(histperfly.frac_linear,1);
histperexp.meanfrac_linear = nan(1,nbins_linear);
histperexp.stdfrac_linear = nan(1,nbins_linear);
for j = 1:nbins_linear,
  [histperexp.meanfrac_linear(j),s] = ...
    weighted_mean_cov(histperfly.frac_linear(j,goodidx)',...
    histperfly.fracframesanalyzed(goodidx)');
  histperexp.stdfrac_linear(j) = sqrt(s);
end
nbins_log = size(histperfly.frac_log,1);
histperexp.meanfrac_log = nan(1,nbins_log);
histperexp.stdfrac_log = nan(1,nbins_log);
for j = 1:size(histperfly.frac_log,1),
  [histperexp.meanfrac_log(j),s] = ...
    weighted_mean_cov(histperfly.frac_log(j,goodidx)',...
    histperfly.fracframesanalyzed(goodidx)');
  histperexp.stdfrac_log(j) = sqrt(s);
end
[histperexp.meanfracless_linear,s] = ...
  weighted_mean_cov(histperfly.fracless_linear(goodidx)',...
  histperfly.fracframesanalyzed(goodidx)');
histperexp.stdfracless_linear = sqrt(s);
[histperexp.meanfracmore_linear,s] = ...
  weighted_mean_cov(histperfly.fracmore_linear(goodidx)',...
  histperfly.fracframesanalyzed(goodidx)');
histperexp.stdfracmore_linear = sqrt(s);
[histperexp.meanfracless_log,s] = ...
  weighted_mean_cov(histperfly.fracless_log(goodidx)',...
  histperfly.fracframesanalyzed(goodidx)');
histperexp.stdfracless_log = sqrt(s);
[histperexp.meanfracmore_log,s] = ...
  weighted_mean_cov(histperfly.fracmore_log(goodidx)',...
  histperfly.fracframesanalyzed(goodidx)');
histperexp.stdfracmore_log = sqrt(s);
histperexp.nflies = nnz(goodidx);