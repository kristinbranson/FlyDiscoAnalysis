function conditions = DecodeConditions(conditionname,conditiondict)

conditions = {};

if strcmpi(conditionname,'any'),
  return;
end
conditionnames = strsplit(conditionname,'_');

conditions = {};
off = 0;
for i = 1:numel(conditionnames),
  conditions(off+1:off+numel(conditiondict.(conditionnames{i}))) = conditiondict.(conditionnames{i});
  off = off + numel(conditiondict.(conditionnames{i}));
end