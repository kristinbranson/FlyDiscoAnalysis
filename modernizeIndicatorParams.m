function result = modernizeIndicatorParams(raw_params)

result = raw_params ;
if ~isfield(result, 'OptogeneticExp') ,
  result.OptogeneticExp = 0 ;
end
