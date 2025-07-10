function [statsperexp] = CombinePerFrameStats2(statsperfly,minZboth,minZfly,field,alldata,alldatanorm,norm_minZboth,norm_minZfly,norm_alldata,norm_alldatanorm)

isnorm = isfield(statsperfly,'norm_mean');

nprctiles = size(statsperfly.prctiles,1);

% which flies we will look at
if isnorm,
  norm_goodidx = ~isnan(statsperfly.norm_Z) & ~isnan(statsperfly.norm_mean) & ...
    ~isinf(statsperfly.norm_mean) & (statsperfly.norm_Z >= norm_minZboth) & (statsperfly.norm_Zfly >= norm_minZfly);
  goodidx = ~isnan(statsperfly.raw_Z) & ~isnan(statsperfly.raw_mean) & ...
    ~isinf(statsperfly.raw_mean) & (statsperfly.raw_Z >= minZboth) & (statsperfly.raw_Zfly >= minZfly);
  both_goodidx = norm_goodidx & goodidx;
else
  goodidx = ~isnan(statsperfly.Z) & ~isnan(statsperfly.mean) & ...
  ~isinf(statsperfly.mean) & (statsperfly.Z >= minZboth) & (statsperfly.Zfly >= minZfly);
end

% total number of frames of data analyzed in this experiment
if isnorm,
  statsperexp.raw_Z = sum(statsperfly.raw_Z(goodidx));
  statsperexp.norm_Z = sum(statsperfly.norm_Z(norm_goodidx));
  % effective number of samples assuming independence -- 1 / (1/m) + (1/n))
  statsperexp.Z = round(1/(1/statsperexp.raw_Z + 1/statsperexp.norm_Z)); 
else
  statsperexp.Z = sum(statsperfly.Z(goodidx));
end

% average number of frames of data analyzed per fly, weighted by weight of
% fly
if isnorm,
  statsperexp.raw_meanZ = weighted_mean(statsperfly.raw_Z(goodidx)',...
    statsperfly.raw_fracframesanalyzed(goodidx)');
  statsperexp.norm_meanZ = weighted_mean(statsperfly.norm_Z(norm_goodidx)',...
    statsperfly.norm_fracframesanalyzed(norm_goodidx)');
  statsperexp.meanZ = round(1/(1/statsperexp.raw_meanZ + 1/statsperexp.norm_meanZ));
else
  statsperexp.meanZ = weighted_mean(statsperfly.Z(goodidx)',...
    statsperfly.fracframesanalyzed(goodidx)');
end

% mean of per-fly means, the total weight for a trajectory is the number of
% frames for which doanalyze_fly is true
if isnorm,
  statsperexp.raw_meanmean = weighted_mean(statsperfly.raw_mean(goodidx)',...
    statsperfly.raw_Zfly(goodidx)');
  statsperexp.norm_meanmean = weighted_mean(statsperfly.norm_mean(norm_goodidx)',...
    statsperfly.norm_Zfly(norm_goodidx)');
  statsperexp.meanmean = statsperexp.raw_meanmean-statsperexp.norm_meanmean;
else
  statsperexp.meanmean = weighted_mean(statsperfly.mean(goodidx)',...
    statsperfly.Zfly(goodidx)');
end

% scatter around the meanmean, but the weighting is slightly different --
% the weight for each trajectory is the *fraction* of frames for which
% doanalyze_fly is true
if isnorm,
  
  [raw_s] = weighted_scatter(statsperfly.raw_mean(goodidx)',...
    statsperexp.raw_meanmean,statsperfly.raw_fracframesanalyzed(goodidx)');
  [norm_s] = weighted_scatter(statsperfly.norm_mean(norm_goodidx)',...
    statsperexp.norm_meanmean,statsperfly.norm_fracframesanalyzed(norm_goodidx)');
  s = weighted_scatter(statsperfly.raw_mean(both_goodidx)'-statsperfly.norm_mean(both_goodidx)',...
    statsperexp.meanmean,statsperfly.fracframesanalyzed(both_goodidx)');
  
  % standard deviation of per-fly means
  statsperexp.raw_stdmean = sqrt(raw_s);
  statsperexp.norm_stdmean = sqrt(norm_s);
  statsperexp.stdmean = sqrt(s);
 
else
  
  [s] = weighted_scatter(statsperfly.mean(goodidx)',...
    statsperexp.meanmean,statsperfly.fracframesanalyzed(goodidx)');
  
  % standard deviation of per-fly means
  statsperexp.stdmean = sqrt(s);
  
end
  
% mean of per-fly standard deviations
if isnorm,
  statsperexp.raw_meanstd = weighted_mean(statsperfly.raw_std(goodidx)',statsperfly.raw_Zfly(goodidx)');
  statsperexp.norm_meanstd = weighted_mean(statsperfly.norm_std(norm_goodidx)',statsperfly.norm_Zfly(norm_goodidx)');
  statsperexp.meanstd = weighted_mean(statsperfly.std(both_goodidx)',statsperfly.Zfly(both_goodidx)');
else
  statsperexp.meanstd = weighted_mean(statsperfly.std(goodidx)',statsperfly.Zfly(goodidx)');
end

% mean and standard deviations of per-fly percentiles

if isnorm,
  
  statsperexp.raw_meanprctiles = nan(1,nprctiles);
  statsperexp.norm_meanprctiles = nan(1,nprctiles);
  statsperexp.meanprctiles = nan(1,nprctiles);
  statsperexp.raw_stdprctiles = nan(1,nprctiles);
  statsperexp.norm_stdprctiles = nan(1,nprctiles);
  statsperexp.stdprctiles = nan(1,nprctiles);
  for j = 1:nprctiles,
    statsperexp.raw_meanprctiles(j) = weighted_mean(statsperfly.raw_prctiles(j,goodidx)',statsperfly.raw_Zfly(goodidx)');
    statsperexp.norm_meanprctiles(j) = weighted_mean(statsperfly.norm_prctiles(j,norm_goodidx)',statsperfly.norm_Zfly(norm_goodidx)');
    statsperexp.meanprctiles(j) = weighted_mean(statsperfly.raw_prctiles(j,both_goodidx)'-statsperfly.norm_prctiles(j,both_goodidx)',statsperfly.Zfly(both_goodidx)');

    s = weighted_scatter(statsperfly.raw_prctiles(j,goodidx)',statsperexp.raw_meanprctiles(j),statsperfly.raw_fracframesanalyzed(goodidx)');
    statsperexp.raw_stdprctiles(j) = sqrt(s);
    s = weighted_scatter(statsperfly.norm_prctiles(j,norm_goodidx)',statsperexp.norm_meanprctiles(j),statsperfly.norm_fracframesanalyzed(norm_goodidx)');
    statsperexp.norm_stdprctiles(j) = sqrt(s);
    s = weighted_scatter(statsperfly.raw_prctiles(j,both_goodidx)'-statsperfly.norm_prctiles(j,both_goodidx)',...
      statsperexp.meanprctiles(j),statsperfly.fracframesanalyzed(both_goodidx)');    
    statsperexp.stdprctiles(j) = sqrt(s);
  end

  
else
  statsperexp.meanprctiles = nan(1,nprctiles);
  statsperexp.stdprctiles = nan(1,nprctiles);
  for j = 1:nprctiles,
    statsperexp.meanprctiles(j) = weighted_mean(statsperfly.prctiles(j,goodidx)',statsperfly.Zfly(goodidx)');
    s = weighted_scatter(statsperfly.prctiles(j,goodidx)',statsperexp.meanprctiles(j),statsperfly.fracframesanalyzed(goodidx)');
    statsperexp.stdprctiles(j) = sqrt(s);
  end
end

if isnorm,
  statsperexp.raw_n = nnz(goodidx);
  statsperexp.norm_n = nnz(norm_goodidx);
  statsperexp.n = nnz(both_goodidx);
else
  statsperexp.n = nnz(goodidx);
end

if ismember(field,{'fractime','boutfreq'}),
  statsperexp.combmean = alldata / alldatanorm;
  statsperexp.combstd = nan;
else
  statsperexp.combmean = nanmean(alldata);
  statsperexp.combstd = nanstd(alldata);
end

if isnorm,

  statsperexp.raw_combmean = statsperexp.combmean;
  statsperexp.raw_combstd = statsperexp.combstd;
  if ismember(field,{'fractime','boutfreq'}),
    statsperexp.norm_combmean = norm_alldata / norm_alldatanorm;
    statsperexp.norm_combstd = nan;
  else
    statsperexp.norm_combmean = nanmean(norm_alldata);
    statsperexp.norm_combstd = nanstd(norm_alldata);
  end
  statsperexp.combmean = statsperexp.raw_combmean - statsperexp.norm_combmean;
  statsperexp.combstd = sqrt(statsperexp.raw_combstd^2 + statsperexp.norm_combstd^2);

end