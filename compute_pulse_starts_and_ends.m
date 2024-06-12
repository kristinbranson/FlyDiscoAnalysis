function [first_sample_index_from_pulse_index, last_sample_index_from_pulse_index] = compute_pulse_starts_and_ends(y)
% y should be a 1xn logical vector

% Extract start+end frames, which are guaranteed to be paired, and in order,
% and start with a pulse-start and end with a pulse-end
delta = diff([0 y 0]) ;  % of length n+1, indexed by edges, including a pre-edge and a post-edge
rising_edge_indices = find(delta == +1) ;
falling_edge_indices = find(delta == -1) ;
first_sample_index_from_pulse_index = rising_edge_indices ;  
  % edge(i) comes just before x(i), so if it's a rising edge then x(i) is the first sample in the pulse
last_sample_index_from_pulse_index = falling_edge_indices - 1 ;  
  % edge(i) comes just after x(i-1), so if it's a falling edge then x(i-1) is the last sample in the pulse
