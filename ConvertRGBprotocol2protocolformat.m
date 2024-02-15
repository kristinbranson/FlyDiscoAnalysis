function [protocol, ledcolor] = ConvertRGBprotocol2protocolformat(RGBprotocol, isLedActive)
% AR 20210311 
% takes in RGB format ledprotocol and prodcues ChR (single
% LED) format protocol
% RGBprotocol = RGB format protocol struct
% countactiveLEDS =  [ 0/1 ; 0/1 ; 0/1 ] 0 = no intensity in any step
% for R/G/B. 1 = nonzero intensity value in at least one step
red = logical([1 0 0]) ;
green = logical([0 1 0]) ;
blue = logical([0 0 1]) ;

protocol = struct('expData',[],'stepNum',[],'intensity',[],'pulseWidthSP',[],'pulsePeriodSP',[],'pulseNum',[], ...
    'offTime',[],'delayTime',[],'iteration',[],'duration',[],'ProtocolData',[],'ProtocolHeader',[]);


protocol.stepNum = RGBprotocol.stepNum;
protocol.duration = RGBprotocol.duration;
protocol.delayTime = RGBprotocol.delayTime;
protocol.expData = RGBprotocol.ProtocolData;
protocol.ProtocolHeader =  RGBprotocol.ProtocolHeader;
protocol.ProtocolData = RGBprotocol.ProtocolData;


if isequal(isLedActive, red)
    protocol.intensity = RGBprotocol.Rintensity;
    protocol.pulseWidthSP = RGBprotocol.RpulseWidth;
    protocol.pulsePeriodSP = RGBprotocol.RpulsePeriod;
    protocol.pulseNum = RGBprotocol.RpulseNum;
    protocol.offTime = RGBprotocol.RoffTime;
    protocol.iteration = RGBprotocol.Riteration;
    ledcolor = 'r';
    
elseif isequal(isLedActive, green)
    protocol.intensity = RGBprotocol.Gintensity;
    protocol.pulseWidthSP = RGBprotocol.GpulseWidth;
    protocol.pulsePeriodSP = RGBprotocol.GpulsePeriod;
    protocol.pulseNum = RGBprotocol.GpulseNum;
    protocol.offTime = RGBprotocol.GoffTime;
    protocol.iteration = RGBprotocol.Giteration;
    ledcolor = 'g';
    
elseif isequal(isLedActive, blue)
    protocol.intensity = RGBprotocol.Bintensity;
    protocol.pulseWidthSP = RGBprotocol.BpulseWidth;
    protocol.pulsePeriodSP = RGBprotocol.BpulsePeriod;
    protocol.pulseNum = RGBprotocol.BpulseNum;
    protocol.offTime = RGBprotocol.BoffTime;
    protocol.iteration = RGBprotocol.Biteration;
    ledcolor = 'b';
end
% new format saves duration in seconds
% data collected before 3/31/21 have duration in seconds, post 3/30 will
% updated to ms. 
% protocol.duration = protocol.duration.*1000;
