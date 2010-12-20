function [stat,varargout] = CollateHistStats(obj,statfun,flies,ns,countsperfly,movie2flies,fly2movie,averaging,statperfly) %#ok<INUSL>

% dimensionality of histogram
nd = ndims(countsperfly) - 1;
% sz = size(countsperfly);
% sz = sz(2:end);

% take the intersection of the specified flies and experiments
allflies_perexp = [movie2flies{ns}];
flies = intersect(flies,allflies_perexp);
ns = unique(fly2movie(flies));
% nflies = length(flies);
% nexpdirs = length(ns);

% is statperfly already computed?
isstatperfly = exist('statperfly','var');

% do we need to compute it?
needstatperfly = nargout > 1 || ismember(lower(averaging),{'allexps_perfly','perexp_perfly'});
if needstatperfly && ~isstatperfly
  Zperfly = sum(countsperfly(:,:),2);
  fracperfly = bsxfun(@rdivide,countsperfly,Zperfly);
  statperfly = statfun(fracperfly);
end

switch lower(averaging),
  
  case 'allexps_allflies',
    
    % treat all frames of data the same: just add up all the counts
    counts = sum(countsperfly(flies,:),1);
    frac = counts / sum(counts(:));
    stat = statfun(frac);
    
  case 'allexps_perfly',
    
    % get per-fly fracs, but treat all flies in all experiments the same
    stat = nanmean(statperfly,1);
    
  case 'perexp_allflies',
    % get per-exp fracs, but treat all frames within an experiment the
    % same
    
    % loop over experiments
    for n = ns(:)',
      % which flies in this experiment
      flies_curr = intersect(flies,movie2flies{n});
      % frac for these flies
      frac_curr = sum(countsperfly(flies_curr,:),1);
      frac_curr = frac_curr / sum(frac_curr(:));
      % compute the statistc
      stat_curr = statfun(frac_curr);
      if n == ns(1),
        % allocate
        stat = zeros(size(stat_curr));
        nexps = stat;
      end
      % only include good data
      isdata = ~isnan(stat_curr);
      stat(isdata) = stat(isdata) + stat_curr(isdata);
      nexps(isdata) = nexps(isdata) + 1;
    end
    % normalize
    stat = stat ./ nexps;
    
  case 'perexp_perfly',
    
    % loop over experiments
    for n = ns(:)',
      
      % which flies in this experiment
      flies_curr = intersect(flies,movie2flies{n});
      
      % average histogram for this experiment
      stat_curr = nanmean(statperfly(flies_curr,:),1);

      if n == ns(1),
        % allocate
        stat = zeros(size(stat_curr));
        nexps = stat;
      end
      
      % only include good data
      isdata = ~isnan(stat_curr);
      stat(isdata) = stat(isdata) + stat_curr(isdata);
      nexps(isdata) = nexps(isdata) + 1;
      
    end
    % normalize
    stat = stat ./ nexps;
    
  otherwise
    
    error('Unknown averaging method %s',averaging);
    
end

if nd > 1,
  stat = shiftdim(stat,1);
end

if nargout > 1,
  varargout = {statperfly};
end