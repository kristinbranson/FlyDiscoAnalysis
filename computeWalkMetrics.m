function [walk_metrics] = computeWalkMetrics(obj,digital_signal,condition_name)


nflies = numel(obj.limbBoutData);

% compute walk metrics for each fly 
for fly = 1:nflies
walkfeaturestruct = compute_WalkFeatures(obj,fly,digital_signal,condition_name);
walk_metrics(fly).walkfeawalkfeaturestruct;




% combine into pairs and all_limbs for each fly 



end


% combine perfly to allflies for perlimb, pairs, and all_limbs