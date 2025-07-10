function histperexp = CombinePerFrameHists2(histperfly,minZboth,minZfly)

% which flies we will look at
goodidx = ~isnan(histperfly.Z) & ~any(isnan(histperfly.frac_linear)) & ...
  ~any(isinf(histperfly.frac_linear)) & ...
   (histperfly.Z >= minZboth) & (histperfly.Zfly >= minZfly);
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
   % mean of per-fly fracs, the total weight for a trajectory is the number of
   % frames for which doanalyze_fly is true
   [histperexp.meanfrac_linear(j)] = ...
     weighted_mean(histperfly.frac_linear(j,goodidx)',...
     histperfly.Zfly(goodidx)');
   % scatter around the meanfrac, but the weighting is slightly different --
   % the weight for each trajectory is the *fraction* of frames for which
   % doanalyze_fly is true
   [s] = weighted_scatter(histperfly.frac_linear(j,goodidx)',...
     histperexp.meanfrac_linear(j),histperfly.fracframesanalyzed(goodidx)');
   histperexp.stdfrac_linear(j) = sqrt(s);
end

nbins_log = size(histperfly.frac_log,1);
histperexp.meanfrac_log = nan(1,nbins_log);
histperexp.stdfrac_log = nan(1,nbins_log);
for j = 1:size(histperfly.frac_log,1),
   % mean of per-fly fracs, the total weight for a trajectory is the number of
   % frames for which doanalyze_fly is true
   [histperexp.meanfrac_log(j)] = ...
     weighted_mean(histperfly.frac_log(j,goodidx)',...
     histperfly.Zfly(goodidx)');
   % scatter around the meanfrac, but the weighting is slightly different --
   % the weight for each trajectory is the *fraction* of frames for which
   % doanalyze_fly is true
   [s] = weighted_scatter(histperfly.frac_log(j,goodidx)',...
     histperexp.meanfrac_log(j),histperfly.fracframesanalyzed(goodidx)');
   histperexp.stdfrac_log(j) = sqrt(s);
end

[histperexp.meanfracless_linear] = ...
  weighted_mean(histperfly.fracless_linear(goodidx)',...
  histperfly.Zfly(goodidx)');
[s] = weighted_scatter(histperfly.fracless_linear(goodidx)',...
  histperexp.meanfracless_linear,histperfly.fracframesanalyzed(goodidx)');
histperexp.stdfracless_linear = sqrt(s);

[histperexp.meanfracmore_linear] = ...
  weighted_mean(histperfly.fracmore_linear(goodidx)',...
  histperfly.Zfly(goodidx)');
[s] = weighted_scatter(histperfly.fracmore_linear(goodidx)',...
  histperexp.meanfracmore_linear,histperfly.fracframesanalyzed(goodidx)');
histperexp.stdfracmore_linear = sqrt(s);

[histperexp.meanfracless_log] = ...
  weighted_mean(histperfly.fracless_log(goodidx)',...
  histperfly.Zfly(goodidx)');
[s] = weighted_scatter(histperfly.fracless_log(goodidx)',...
  histperexp.meanfracless_log,histperfly.fracframesanalyzed(goodidx)');
histperexp.stdfracless_log = sqrt(s);

[histperexp.meanfracmore_log] = ...
  weighted_mean(histperfly.fracmore_log(goodidx)',...
  histperfly.Zfly(goodidx)');
[s] = weighted_scatter(histperfly.fracmore_log(goodidx)',...
  histperexp.meanfracmore_log,histperfly.fracframesanalyzed(goodidx)');
histperexp.stdfracmore_log = sqrt(s);

histperexp.n = nnz(goodidx);