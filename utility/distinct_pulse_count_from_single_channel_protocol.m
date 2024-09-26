function result = distinct_pulse_count_from_single_channel_protocol(protocol, n_pad)
% Count the number of distinct pulses in the protocol.
% "Distinct" means there has to be a gap between them.
y_raw = signal_from_single_channel_protocol(protocol) ;
y = compute_pulse_envelope(y_raw, n_pad) ;
y_binary = (y>0) ;
dy_binary = diff(y_binary) ;
is_rising_edge = (dy_binary>0) ;
result = sum(is_rising_edge) ;
end
