function conditions = DecodeConditions(conditionname,conditiondict)

conditions = {};

if strcmpi(conditionname,'any'),
  return;
end
conditionnames = strsplit(conditionname,'_');

conditions = cell(1,2*numel(conditionnames));
for i = 1:numel(conditionnames),
  conditions(2*i-1:2*i) = conditiondict.(conditionnames{i});
end