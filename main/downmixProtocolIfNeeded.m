function result = downmixProtocolIfNeeded(raw_protocol, indicator_params)
% Examine the metadata and the contents of the LED protocol file to determine
% the protocol.  A protocol from an RGB experiment has to be massaged a bit.
% Other protocols are passed through unaltered.

% Deal with arguments
if ~exist('indicator_params', 'var') || isempty(indicator_params) ,
  indicator_params = [] ;
end

if isfield(raw_protocol, 'intensity') ,
  % This looks like a traditional one-LED protocol
  result = raw_protocol ;
else
  has_explicit_masks =  ~isempty(indicator_params) && isfield(indicator_params, 'protocol_prefix_from_mask_index') ;
  if has_explicit_masks ,
    protocol_prefix_from_mask_index = indicator_params.protocol_prefix_from_mask_index ;
    primary_channel_prefix = protocol_prefix_from_mask_index{1} ;
  else
    % If protocol is missing the itensity field, it is presumably an RGB protocol.  
    % But with no explicit maks, just have to guess about what the primary
    % indicator is.
    primary_channel_prefix = guessPrimaryChannelProtocolPrefixFromRawProtocol(raw_protocol) ;
  end
  % call function that transforms new protocol to old protocol
  result = single_channel_protocol_from_multi_channel_protocol(raw_protocol, primary_channel_prefix) ;
end
