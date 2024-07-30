function result = downmix_indicatorLED(indicatorLED)
% Take the indicatorLED struct, and make it look like indicatorLED looked
% before we added proper multi-LED supprt.  We always use the 'primary'
% indicator, because that's the one that should carry the 'primary' stimulus
% information.

result = indicatorLED ;
if iscell(indicatorLED.startframe) ,
  if isempty(indicatorLED.startframe) ,
    error('indicatorLED.startframe is empty.  Likely there was a problem with indicator detection.') ;
  else
    % Look at the "primary" indicator
    result.startframe = indicatorLED.startframe{1} ;
    result.endframe = indicatorLED.endframe{1} ;
    result.starttimes = indicatorLED.starttimes{1} ;
    result.endtimes = indicatorLED.endtimes{1} ;
  end
end

end  % function
