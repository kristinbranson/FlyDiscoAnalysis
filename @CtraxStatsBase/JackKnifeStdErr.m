function [jackknife_stat_stderr,jackknife_stat_mean] = ...
  JackKnifeStdErr(statfun,nflies,nexpdirs,jackknife)

jackknife_stat_mean = [];
jackknife_stat_stderr = [];

if strcmpi(jackknife,'none'),
  return;
end

switch jackknife,

  case 'perfly',
    for fly = 1:nflies,
      
      % use all but the current fly
      fliescurr = setdiff(1:nflies,fly);
      
      stat_curr = statfun(fliescurr,1:nexpdirs);
      
      % allocate
      if fly == 1,
        sz = size(stat_curr);
        jackknife_stat_mean = zeros(sz);
        jackknife_stat_stderr = zeros(sz);
        nsets = zeros(sz);
      end
      
      % make sure not nan
      isdata = ~isnan(stat_curr);
      nsets(isdata) = nsets(isdata) + 1;
      
      % update mean, std estimates
      jackknife_stat_mean(isdata) = jackknife_stat_mean(isdata) + stat_curr(isdata); %#ok<*AGROW>
      jackknife_stat_stderr(isdata) = jackknife_stat_stderr(isdata) + stat_curr(isdata).^2;
    end
    
  case 'perexp',
  
    for n = 1:nexpdirs,
      
      % use all but the current exp
      nscurr = setdiff(1:nexpdirs,n);

      stat_curr = statfun(1:nflies,nscurr);
      
      % allocate
      if n == 1,
        sz = size(stat_curr);
        jackknife_stat_mean = zeros(sz);
        jackknife_stat_stderr = zeros(sz);
        nsets = zeros(sz);
      end

      % make sure not nan
      isdata = ~isnan(stat_curr);
      nsets(isdata) = nsets(isdata) + 1;
      
      % update mean, stderr estimates
      jackknife_stat_mean(isdata) = jackknife_stat_mean(isdata) + stat_curr(isdata);
      jackknife_stat_stderr(isdata) = jackknife_stat_stderr(isdata) + stat_curr(isdata).^2;
      
    end
    
  otherwise
    
    error('Unknown jackknife bagging unit %s',jackknife);
    
end
% normalize
jackknife_stat_mean = jackknife_stat_mean ./ nsets;
jackknife_stat_stderr = jackknife_stat_stderr ./ nsets - jackknife_stat_mean.^2;
jackknife_stat_stderr = sqrt(jackknife_stat_stderr .* (nsets-1)./nsets);
