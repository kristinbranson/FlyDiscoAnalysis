function [result, t] = computeRgbStimulusProtocol(protocol)

% Extract out the protocol for each color channel
red_protocol = single_channel_protocol_from_multi_channel_protocol(protocol, 'R') ;
green_protocol = single_channel_protocol_from_multi_channel_protocol(protocol, 'G') ;
blue_protocol = single_channel_protocol_from_multi_channel_protocol(protocol, 'B') ;

% Compute the signal for each color channel
dt = 0.001 ;  % s
[Yr,t] = signal_from_single_channel_protocol(red_protocol, dt) ;
[Yg,~] = signal_from_single_channel_protocol(green_protocol, dt) ;
[Yb,~] = signal_from_single_channel_protocol(blue_protocol, dt) ;

% Package them up for return
result = [Yr Yg Yb] ;
