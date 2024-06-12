function LoadIndicatorDataFromFile(obj,n)

try
  indicatordatafile = fullfile(obj.expdirs{n},obj.dataloc_params.indicatordatafilestr);
  if ~exist(indicatordatafile,'file'),
    return;
  end
  id = load(indicatordatafile);
catch ME
  warning(getReport(ME));
  return;
end
indicatordigital = id.indicatorLED.indicatordigital ;
indicatorLED = convertIndicatorData(indicatordigital) ;
obj.indicatorLED{n} = indicatorLED;

end  % function



function result = convertIndicatorData(indicatordigital)
% Extracts the pulse starts/ends from the indicatordigital signal and alose
% the "after-pulse" starts/ends.  And puts it in a struct of the type that the
% rest of the code expects.

% Just easier to recompute this stuff here from the digital signal
% Also, pre-mid-2023 versions of the LED indicator stage have an off-by-one bug,
% so better just to do it here to be sure.
y = logical(indicatordigital) ;
[first_sample_index_from_pulse_index, last_sample_index_from_pulse_index] = compute_pulse_starts_and_ends(y) ;
n = numel(y) ;

% Assemble the result struct.
%
% A "pulse" is a continuous span of high samples, with a low sample just
% before and just after.  (Where samples outside the interval [1,n] are
% considered low.)
%
% An "after-pulse" is the continuous span of low samples in-between two
% pulses.
result = struct() ;
result.starton = first_sample_index_from_pulse_index ;  % frame index of the first high frame of each pulse 
result.endon = last_sample_index_from_pulse_index ;  % frame index of the last high frame of each pulse
result.startoff = last_sample_index_from_pulse_index+1 ;  % frame index of the first low frame of each after-pulse 
result.endoff = [first_sample_index_from_pulse_index(2:end)-1 n] ;  % frame index of the last low frame of each after-pulse

% % Debug code
% pulse_count = numel(result.starton) ;
% figure('color', 'w') ;
% plot(y,'b.') ;
% ylim([-0.05 1.05]) ;
% hold on ;
% plot(result.starton, ones(1, pulse_count), 'og', 'MarkerSize', 9) ;
% plot(result.endon, ones(1, pulse_count), 'or', 'MarkerSize', 9) ;
% plot(result.startoff, zeros(1, pulse_count), 'ob', 'MarkerSize', 9) ;
% plot(result.endoff, zeros(1, pulse_count), 'om', 'MarkerSize', 9) ;
% hold off; 

% Quick sanity-check
assert(numel(result.starton) == numel(result.endon)) ;
assert(numel(result.startoff) == numel(result.endoff)) ;

end  % function
