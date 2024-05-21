function result = distinct_pulse_count_from_single_channel_protocol(protocol)
% Count the number of distinct pulses in the protocol.
% "Distinct" means there has to be a gap between them.
y = signal_from_single_channel_protocol(protocol) ;
dy = diff(y) ;
is_rising_edge = (dy>0) ;
result = sum(is_rising_edge) ;
end