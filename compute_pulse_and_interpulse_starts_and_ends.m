function [first_tick_index_from_pulse_index, ...
          last_tick_index_from_pulse_index, ...
          first_tick_index_from_interpulse_index, ...
          last_tick_index_from_interpulse_index] = ...
  compute_pulse_and_interpulse_starts_and_ends(y_from_tick_index)
% Compute the start and end of each pulse and "interpulse" in the digital
% signal y_from_tick_index, which should be logical row vector.  A "pulse" is
% a contiguous run of highs, an "interpulse" a contiguous run of lows. Note
% that a run of lows at the start or end of the signal is considered an
% interpulse.  Also note that the return values are such that:
%
%   all(y_from_tick_index(first_tick_index_from_pulse_index))
%   all(y_from_tick_index(last_tick_from_pulse_index))
%   all(~y_from_tick_index(first_tick_index_from_interpulse_index))
%   all(~y_from_tick_index(last_tick_from_interpulse_index))

% Get the size of the input  
tick_count = numel(y_from_tick_index) ;

% Handle an empty input as a special case
if tick_count == 0 ,
  first_tick_index_from_pulse_index = zeros(1,0) ;
  last_tick_index_from_pulse_index = zeros(1,0) ;
  first_tick_index_from_interpulse_index = zeros(1,0) ;
  last_tick_index_from_interpulse_index = zeros(1,0) ;
  return
end

% If we get here, tick_count >= 1

% Pre-allocate these for speed
first_tick_index_from_pulse_index = zeros(1,tick_count) ;
last_tick_index_from_pulse_index = zeros(1,tick_count) ;
first_tick_index_from_interpulse_index = zeros(1,tick_count) ;
last_tick_index_from_interpulse_index = zeros(1,tick_count) ;

% No pulses or interpulses yet
pulse_count = 0 ;
interpulse_count = 0 ;

% Deal with the first element
tick_index = 1 ;
y = y_from_tick_index(tick_index) ;
if y ,
  % This is the start of the first pulse
  pulse_count = pulse_count + 1 ;
  pulse_index = pulse_count ;
  first_tick_index_from_pulse_index(pulse_index) = tick_index ;
else
  % This is the start of the first interpulse
  interpulse_count = interpulse_count + 1 ;
  interpulse_index = interpulse_count ;
  first_tick_index_from_interpulse_index(interpulse_index) = tick_index ;
end
% Initialize y_last
y_last = y ;

% Iterate through the rest of the elements (JIT accelerator should make this
% fast)
for tick_index = 2 : tick_count ,
  y = y_from_tick_index(tick_index) ;
  if y ~= y_last ,
    if y ,
      % This is a rising edge
      % Finish the interpulse
      last_tick_index_from_interpulse_index(interpulse_index) = tick_index-1 ;
      % Deal with the new pulse
      pulse_count = pulse_count + 1 ;
      pulse_index = pulse_count ;
      first_tick_index_from_pulse_index(pulse_index) = tick_index ;
    else
      % This is a falling edge
      % Finish the pulse
      last_tick_index_from_pulse_index(pulse_index) = tick_index-1 ;
      % Deal with the new interpulse
      interpulse_count = interpulse_count + 1 ;
      interpulse_index = interpulse_count ;
      first_tick_index_from_interpulse_index(interpulse_index) = tick_index ;
    end
    % Update y_last
    y_last = y ;    
  end
end

% Now that we've reached the end, need to finish off the current
% pulse/interpulse
if y_last ,
  % finish the pulse
  last_tick_index_from_pulse_index(pulse_index) = tick_count ;  
else
  % finish the interpulse
  last_tick_index_from_interpulse_index(interpulse_index) = tick_count ;
end

% Right-size the preallocated arrays before returning
first_tick_index_from_pulse_index = first_tick_index_from_pulse_index(1:pulse_count) ;
last_tick_index_from_pulse_index = last_tick_index_from_pulse_index(1:pulse_count) ;
first_tick_index_from_interpulse_index = first_tick_index_from_interpulse_index(1:interpulse_count) ;
last_tick_index_from_interpulse_index = last_tick_index_from_interpulse_index(1:interpulse_count) ;
