function [are_timestamps_reliable, fallback_dt] = getFallbackDtIfNeeded(registration_params, timestamps, trx)
% Compute the median interval, taking into account the registration_params,
% timestamps, and trx.  
% This is a pure function.

are_timestamps_reliable = registration_params.usemediandt ;  % This field is maybe suboptimally named?
if are_timestamps_reliable ,
  fallback_dt = [] ;  % this should not be used if timestamps are reliable
else
  dt = diff(timestamps);
  dt(isnan(dt)) = [];
  if isempty(dt),
    if ~isempty(trx) ,
      fallback_dt = 1/trx(1).fps;
    else
      error('Cannot compute median interval because there are no non-nan dt''s from timestamps and trx is empty') ;
    end
  else
    fallback_dt = median(dt);
  end
end

end
