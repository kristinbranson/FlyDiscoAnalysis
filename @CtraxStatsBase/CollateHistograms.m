function [frac,varargout] = CollateHistograms(flies,ns,countsperfly,movie2flies,fly2movie,averaging,fracperfly)

% dimensionality of histogram
nd = ndims(countsperfly) - 1;
sz = size(countsperfly);
sz = sz(2:end);

% take the intersection of the specified flies and experiments
allflies_perexp = [movie2flies{ns}];
flies = intersect(flies,allflies_perexp);
ns = unique(fly2movie(flies));
% nflies = length(flies);
% nexpdirs = length(ns);

% is fracperfly already computed?
isfracperfly = exist('fracperfly','var');

% do we need to compute it?
needfracperfly = nargout > 1 || ismember(lower(averaging),{'allexps_perfly','perexp_perfly'});
if needfracperfly && ~isfracperfly
  Zperfly = sum(countsperfly(:,:),2);
  fracperfly = bsxfun(@rdivide,countsperfly,Zperfly);
end

switch lower(averaging),
  
  case 'allexps_allflies',
    
    % treat all frames of data the same: just add up all the counts
    counts = sum(countsperfly(flies,:),1);
    frac = counts / sum(counts(:));
    
  case 'allexps_perfly',
    
    % get per-fly fracs, but treat all flies in all experiments the same
    frac = nanmean(fracperfly(flies,:),1);
    
  case 'perexp_allflies',
    % get per-exp fracs, but treat all frames within an experiment the
    % same
    
    % loop over experiments
    frac = zeros([1,sz]);
    nexps = zeros([1,sz]);
    for n = ns(:)',
      % which flies in this experiment
      flies_curr = intersect(flies,movie2flies{n});
      % frac for these flies
      frac_curr = sum(countsperfly(flies_curr,:),1);
      frac_curr = frac_curr / sum(frac_curr(:));
      % only include good data
      isdata = ~isnan(frac_curr);
      frac(isdata) = frac(isdata) + frac_curr(isdata);
      nexps(isdata) = nexps(isdata) + 1;
    end
    % normalize
    frac = frac ./ nexps;
    
  case 'perexp_perfly',
    
    % loop over experiments
    frac = zeros([1,sz]);
    nexps = zeros([1,sz]);
    for n = ns(:)',
      
      % which flies in this experiment
      flies_curr = intersect(flies,movie2flies{n});
      
      % average histogram for this experiment
      frac_curr = nanmean(fracperfly(flies_curr,:),1);
      
      % only include good data
      isdata = ~isnan(frac_curr);
      frac(isdata) = frac(isdata) + frac_curr(isdata);
      nexps(isdata) = nexps(isdata) + 1;
      
    end
    % normalize
    frac = frac ./ nexps;
    
  otherwise
    
    error('Unknown averaging method %s',averaging);
    
end

if nd > 1,
  frac = reshape(frac,sz);
end

if nargout > 1,
  varargout = {fracperfly};
end