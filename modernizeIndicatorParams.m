function result = modernizeIndicatorParams(raw_params)

result = raw_params ;
if ~isfield(raw_params, 'OptogeneticExp') ,
  result.OptogeneticExp = 0 ;
end
% If there's only one channel, protocol_prefix_from_mask_index will come
% through as an old-style string, but we want a cell array of old-style
% strings.
if isfield(raw_params, 'protocol_prefix_from_mask_index') && ischar(raw_params.protocol_prefix_from_mask_index) ,
  result.protocol_prefix_from_mask_index = { raw_params.protocol_prefix_from_mask_index } ;
end
