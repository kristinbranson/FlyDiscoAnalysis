function meddt = medianIntervalFromTimestamps(timestamps)
% Compute the median interval from an array of timestamps.
% This is a pure function.

dt = diff(timestamps);
dt(isnan(dt)) = [];
if isempty(dt),
  meddt = 1/trx(1).fps;
else
  meddt = median(dt);
end

end
