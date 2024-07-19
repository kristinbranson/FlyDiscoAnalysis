function result = do_run_from_analysis_parameters(protocol)

% Convert individual fields in protocol with fields so<stage_name>
% to fields named <stage_name>, and with values coerced to one of 'on', 'off',
% 'force'.

stage_names = FlyDiscoStageNames() ;
stage_count = numel(stage_names) ;
result = struct() ;
for stage_index = 1 : stage_count ,
  stage_name = stage_names{stage_index} ;
  protocol_field_name = sprintf('do%s', stage_name) ;
  raw_value = protocol.(protocol_field_name) ;
  result.(stage_name) = coerce_to_on_off_force(raw_value) ;
end

end
