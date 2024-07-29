function result = single_channel_protocol_from_multi_channel_protocol(input_protocol, primary_channel_prefix)
% Convert a multi-channel protocol to a single-channel protocol, given the
% 'prefix' of the primary channel.  For instance, if the primary channel has
% field names like 'Rintensity', 'RpulseWidth', etc, then the
% primary_channel_prefix would be 'R'.

result = struct();
result.stepNum = input_protocol.stepNum ;
result.duration = input_protocol.duration ;
result.delayTime = input_protocol.delayTime ;
result.intensity = input_protocol.(strcat(primary_channel_prefix, 'intensity')) ;
result.pulseWidthSP = input_protocol.(strcat(primary_channel_prefix, 'pulseWidth')) ;
result.pulsePeriodSP = input_protocol.(strcat(primary_channel_prefix, 'pulsePeriod')) ;
result.pulseNum = input_protocol.(strcat(primary_channel_prefix, 'pulseNum')) ;
result.offTime = input_protocol.(strcat(primary_channel_prefix, 'offTime')) ;
result.iteration = input_protocol.(strcat(primary_channel_prefix, 'iteration')) ;
result.expData = input_protocol.ProtocolData ;
result.ProtocolHeader =  input_protocol.ProtocolHeader ;
result.ProtocolData = input_protocol.ProtocolData ;
