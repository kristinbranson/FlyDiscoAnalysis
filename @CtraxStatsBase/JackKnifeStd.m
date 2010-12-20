function [jackknife_stat_std,jackknife_stat_mean] = ...
  JackKnifeStd(statfun,nflies,nexpdirs,jackknife)

jackknife_stat_mean = [];
jackknife_stat_std = [];

if strcmpi(jackknife,'none'),
  return;
end

switch jackknife,

  case 'perfly',
    for fly = 1:nflies,
      % use only the current fly
      stat_curr = statfun(fly,1:nexpdirs);
      
      % allocate
      if fly == 1,
        sz = size(stat_curr);
        jackknife_stat_mean = zeros(sz);
        jackknife_stat_std = zeros(sz);
        nsets = zeros(sz);
      end
      
      % make sure not nan
      isdata = ~isnan(stat_curr);
      nsets(isdata) = nsets(isdata) + 1;
      
      % update mean, std estimates
      jackknife_stat_mean(isdata) = jackknife_stat_mean(isdata) + stat_curr(isdata); %#ok<*AGROW>
      jackknife_stat_std(isdata) = jackknife_stat_std(isdata) + stat_curr(isdata).^2;
    end
    
  case 'perexp',
  
    for n = 1:nexpdirs,
      
      % use only the current exp
      stat_curr = statfun(1:nflies,n);
      
      % allocate
      if n == 1,
        sz = size(stat_curr);
        jackknife_stat_mean = zeros(sz);
        jackknife_stat_std = zeros(sz);
        nsets = zeros(sz);
      end

      % make sure not nan
      isdata = ~isnan(stat_curr);
      nsets(isdata) = nsets(isdata) + 1;
      
      % update mean, stderr estimates
      jackknife_stat_mean(isdata) = jackknife_stat_mean(isdata) + stat_curr(isdata);
      jackknife_stat_std(isdata) = jackknife_stat_std(isdata) + stat_curr(isdata).^2;
      
    end
    
  otherwise
    
    error('Unknown jackknife bagging unit %s',jackknife);
    
end
% normalize
jackknife_stat_mean = jackknife_stat_mean ./ nsets;
jackknife_stat_std = sqrt(jackknife_stat_std ./ nsets - jackknife_stat_mean.^2);

