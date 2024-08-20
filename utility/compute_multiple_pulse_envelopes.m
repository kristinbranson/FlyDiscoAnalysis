function [y, first_sample_index_from_pulse_index, last_sample_index_from_pulse_index] = compute_multiple_pulse_envelopes(x, n_pad)
% Like compute_pulse_envelope(), but for n_signals x n_t inputs.
% See docs for compute_pulse_envelope().  On return, y is the same shape as x,
% first_sample_index_from_pulse_index is a cell array of shape n_signals x 1,
% each element holding a 1 x pulse_count double array, and
% last_sample_index_from_pulse_index is similar to
% first_sample_index_from_pulse_index.

signal_count = size(x,1) ;
tick_count = size(x,2) ;
y = false(signal_count, tick_count) ;
first_sample_index_from_pulse_index = cell(signal_count, 1) ;
last_sample_index_from_pulse_index = cell(signal_count, 1) ;
for j = 1 : signal_count ,
  [y(j,:), first_sample_index_from_pulse_index{j}, last_sample_index_from_pulse_index{j}] = ...
    compute_pulse_envelope(x(j,:), n_pad) ;
end
