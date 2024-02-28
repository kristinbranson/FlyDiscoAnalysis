function result = downmixProtocolIfNeeded(rawProtocol)
% Examine the metadata and the contents of the LED protocol file to determine
% the protocol.  A protocol from an RGB experiment has be massaged a bit.
% Other protocols are passed through unaltered.

if isfield(rawProtocol, 'intensity') ,
  % This looks like a traditional one-LED protocol
  result = rawProtocol ;
else
  % If protocol is missing the itensity field, it is presumably an RGB protocol.
  % But check to make sure.
  if ~isfield(rawProtocol, 'Rintensity') || ~isfield(rawProtocol, 'Gintensity') || ~isfield(rawProtocol, 'Bintensity') ,
    error('Protocol lacks intensity field, but is also missing at least one of Rintensity, Gintensity, Bintensity fields.') ;
  end
  % Test if RGBprotocol has only one active color
  isActiveFromLedIndex = [ any(rawProtocol.Rintensity) any(rawProtocol.Gintensity) any(rawProtocol.Bintensity) ] ;
  % Currently we only use the first active LED.  Warn if there is more than one.
  activeLedIndices = find(isActiveFromLedIndex) ;
  if isempty(activeLedIndices) ,
    error('RGB protocol has zero active LEDs')
  elseif length(activeLedIndices)>1 ,
    warning('RGB protocol has more than one active LED---only using the first one')    
  end
  willUseLedIndex = activeLedIndices(1) ;
  willUseFromLedIndex = ([1 2 3]==willUseLedIndex) ;
  % call function that transforms new protocol to old protocol
  result = ConvertRGBprotocol2protocolformat(rawProtocol, willUseFromLedIndex) ;
end
