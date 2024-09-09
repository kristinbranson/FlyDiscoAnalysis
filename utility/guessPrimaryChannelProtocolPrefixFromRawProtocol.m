function result = guessPrimaryChannelProtocolPrefixFromRawProtocol(rawProtocol)
  % In an RGB protocol, there are typically protocols for multiple channels.
  % For instance, the available intensity fields might be named Rintensity,
  % Gintensity, and Bintensity.  Without more information, it's hard to say
  % which of these corresponds to the primary stimulus channel.  This function
  % assumes the primary channel is one of R, G, and B, determines which of these
  % channels are 'active', and then returns the 'prefix' (one of 'R', 'G', and
  % 'B') of the first active channel.  A channel is active if it has at least
  % one train with nonzero intensity.

  % Detetmine which of the channels are active
  prefixFromChannelIndex = { 'R', 'G', 'B' } ;  
  prefixCount = numel(prefixFromChannelIndex) ;
  isActiveFromChannelIndex = false(1,prefixCount) ;
  for prefixIndex = 1 : prefixCount ,
    prefix = prefixFromChannelIndex{prefixIndex} ;
    fieldName = horzcat(prefix, 'intensity') ;
    if isfield(rawProtocol, fieldName) ,
      isActiveFromChannelIndex(prefixIndex) = any(rawProtocol.(fieldName)) ;
    end
  end

  % Determine the first active channel.  Warn if there is more than one.
  activeChannelIndices = find(isActiveFromChannelIndex) ;
  if isempty(activeChannelIndices) ,
    error('For experiment with no explicit masks, RGB protocol has zero active LEDs')
  elseif isscalar(activeChannelIndices) ,
    primaryChannelIndex = activeChannelIndices ;       
  else 
    primaryChannelIndex = activeChannelIndices(1) ;
    warning('For experiment with no explicit masks, RGB protocol has more than one active LED---assuming the first active LED (LED %d) is primary', ...
            primaryChannelIndex)
  end
  result = prefixFromChannelIndex{primaryChannelIndex} ;
end
