function result = downmixProtocolIfNeeded(metadata, rawProtocol)
% Examine the metadata and the contents of the LED protocol file to determine
% the protocol.  A protocol from an RGB experiment has be massaged a bit.
% Other protocols are passed through unaltered.

if isExperimentRGB(metadata) && isfield(rawProtocol,'Rintensity')
  % test if RGBprotocol has only one active color
  isActiveFromLedIndex = [ any(rawProtocol.Rintensity) any(rawProtocol.Gintensity) any(rawProtocol.Bintensity) ] ;
  % check that there is 1 and only 1 color LED used in protocol
  activeLedCount = sum(double(isActiveFromLedIndex)) ;
  if activeLedCount == 0 ,
    error('ChR = 1 for LED protcol with no active LEDs')
  elseif activeLedCount == 1 ,
    % do nothing, this is what we want
  else
    error('More than one active LED color in protocol. Not currently supported')
  end
  % call function that transforms new protocol to old protocol
  result = ConvertRGBprotocol2protocolformat(rawProtocol, isActiveFromLedIndex) ;
else
  result = rawProtocol ;
end
