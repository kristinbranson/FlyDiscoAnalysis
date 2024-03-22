function result = downmixProtocolIfNeeded(rawProtocol)
% Examine the metadata and the contents of the LED protocol file to determine
% the protocol.  A protocol from an RGB experiment has be massaged a bit.
% Other protocols are passed through unaltered.

if isfield(rawProtocol, 'intensity') ,
  % This looks like a traditional one-LED protocol
  result = rawProtocol ;
else
  % If protocol is missing the itensity field, it is presumably an RGB protocol.  
  isActiveFromLedIndex = false(1,3) ;
  if isfield(rawProtocol, 'Rintensity') ,
    isActiveFromLedIndex(1) = any(rawProtocol.Rintensity) ;
  end
  if isfield(rawProtocol, 'Gintensity') ,
    isActiveFromLedIndex(2) = any(rawProtocol.Gintensity) ;
  end
  if isfield(rawProtocol, 'Bintensity') ,
    isActiveFromLedIndex(3) = any(rawProtocol.Bintensity) ;
  end
  % Currently we only use the first active LED.  Warn if there is more than one.
  % TODO: This code should be removed once we're doing proper handling of RGB
  % protocols.  -- ALT, 2024-02-29
  activeLedIndices = find(isActiveFromLedIndex) ;
  if isempty(activeLedIndices) ,
    error('RGB protocol has zero active LEDs')
  elseif isscalar(activeLedIndices) ,
    willUseLedIndex = activeLedIndices ;       
  else 
    willUseLedIndex = activeLedIndices(1) ;
    warning('RGB protocol has more than one active LED---only using the first active LED (LED %d)', willUseLedIndex)
  end
  willUseFromLedIndex = ([1 2 3]==willUseLedIndex) ;
  % call function that transforms new protocol to old protocol
  result = ConvertRGBprotocol2protocolformat(rawProtocol, willUseFromLedIndex) ;
end
