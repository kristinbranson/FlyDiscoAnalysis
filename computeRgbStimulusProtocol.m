function [result, t] = computeRgbStimulusProtocol(protocol)

% Extract out the protocol for each color channel
red_protocol = ConvertRGBprotocol2protocolformat(protocol, [1 0 0]) ;
green_protocol = ConvertRGBprotocol2protocolformat(protocol, [0 1 0]) ;
blue_protocol = ConvertRGBprotocol2protocolformat(protocol, [0 0 1]) ;

% Compute the signal for each color channel
dt = 0.001 ;  % s
[Yr,t] = signal_from_single_channel_protocol(red_protocol, dt) ;
[Yg,~] = signal_from_single_channel_protocol(green_protocol, dt) ;
[Yb,~] = signal_from_single_channel_protocol(blue_protocol, dt) ;

% Package them up for return
result = [Yr Yg Yb] ;
