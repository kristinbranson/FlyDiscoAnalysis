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
indicatorLed = LoadIndicatorDataFromFileCore(id) ;
obj.indicatorLED{n} = indicatorLed;

end  % function





function result = LoadIndicatorDataFromFileCore(id)
% Extracts the pulse starts/ends from the indicatordigital signal and alose
% the "after-pulse" starts/ends.  And puts it in a struct of the type that the
% rest of the code expects.
indicatordigital = id.indicatorLED.indicatordigital(1,:);


% Just easier to recompute this stuff here from the digital signal
% Also, pre-mid-2023 versions of the LED indicator stage have an off-by-one bug,
% so better just to do it here to be sure.
y = logical(indicatordigital) ;  % Just to make sure it's logical
[first_tick_index_from_pulse_index, last_tick_index_from_pulse_index, ...
 first_tick_index_from_interpulse_index, last_tick_index_from_interpulse_index] = ...
  compute_pulse_and_interpulse_starts_and_ends(y) ;

result.starton = first_tick_index_from_pulse_index ;  % frame index of the first high frame of each pulse 
result.endon = last_tick_index_from_pulse_index ;  % frame index of the last high frame of each pulse
result.startoff = first_tick_index_from_interpulse_index ;  % frame index of the first low frame of each interpulse 
result.endoff = last_tick_index_from_interpulse_index ;  % frame index of the last low frame of each interpulse

result.indicatordigital = indicatordigital; % return the full movie indicator digital signal AR

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
