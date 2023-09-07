function [y, first_sample_index_from_pulse_index, last_sample_index_from_pulse_index] = compute_pulse_envelope(x, n_pad)
% y is a boolean signal of the same length as x, indicating when x(i) is
% 'robustly' high, bridging over small gaps in the high state.
% Gaps smaller than n_pad elements are ignored.
%
% first_sample_index_from_pulse_index(pulse_index) is the index of the first sample of x
% within pulse pulse_index.  It is guaranteed to be on [1, n], where n is the
% length of x.
%
% last_sample_index_from_pulse_index(pulse_index) is the index of the last sample of x
% within pulse pulse_index.  It is guaranteed to be on [1, n], where n is the
% length of x.
%
% Things guaranteed to be true on exit:
%   1. length(first_sample_index_from_pulse_index) ==
%      length(last_sample_index_from_pulse_index)
%   2. If i<j, first_sample_index_from_pulse_index(i) <
%      first_sample_index_from_pulse_index(j)
%   3. If i<j, last_sample_index_from_pulse_index(i) <
%      last_sample_index_from_pulse_index(j)
%   4. first_sample_index_from_pulse_index(i) <
%      last_sample_index_from_pulse_index(i), for all i
%
%   That is, the edges start with a rising edge, alternate between falling/rising after 
%   that, and end with a falling pulse. 

% % For each sample, compute the sum of the next n_pad samples
% next_kernel = [ones(1, n_pad) 0 zeros(1, n_pad)] ;
% y_next = conv(y_raw, next_kernel, 'same') 
% 
% % For each sample, compute the sum of the previous n_pad samples
% previous_kernel = [zeros(1, n_pad) 0 ones(1, n_pad)] ;
% y_previous = conv(y_raw, previous_kernel, 'same') 

% For each sample, compute the sum of the previous, next n_pad samples
% This code is faster than the commented-out version above.
kernel = ones(1, n_pad) ;
y_near = conv(x, kernel, 'full') ;
y_next = logical([y_near(n_pad+1:end) 0]) ;
y_previous = logical([0 y_near(1:end-n_pad)]) ;

% Go through the elements, bridging over too-small gaps
% Matlab accelerator should make this loop fast
n = length(x) ;
y = false(1, n) ;
state = false ;
for i = 1:n ,
  value = x(i) ;
  if state ,
    % state is true
    if value ,
      % state is true, value is true
      % No change
    else
      % state is true, value is false
      % Possible falling edge
      if y_next(i) ,
        % There are highs in the near future, therefore:
        % Ignore the falling edge
      else
        % There are all lows in the near future, therefore:
        % Accept the falling edge
        state = false ;
      end
    end
  else
    % state is false
    if value ,
      % state is false, value is true
      % Possible rising edge
      if y_previous(i) ,
        % There have been highs in the recent past, therefore:
        % Ignore the rising edge
      else
        % There are all lows in the recent past, therefore:
        % Accept the rising edge
        state = true ;
      end
    else
      % state is false, value is false
      % No change
    end    
  end
  y(i) = state ;
end

% Extract start+end frames, which are guaranteed to be 
% paired, and in order, and start with start and end with an end
delta = diff([0 y 0]) ;  % of length n+1, indexed by edges, including a pre-edge and a post-edge
rising_edge_indices = find(delta == +1) ;
falling_edge_indices = find(delta == -1) ;
first_sample_index_from_pulse_index = rising_edge_indices ;  
  % edge(i) comes just before x(i), so if it's a rising edge then x(i) is the first sample in the pulse
last_sample_index_from_pulse_index = falling_edge_indices - 1 ;  
  % edge(i) comes just after x(i-1), so if it's a falling edge then x(i-1) is the last sample in the pulse
