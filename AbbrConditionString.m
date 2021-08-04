function conditions = AbbrConditionString(frameconditions,flyconditions)

conditions = cell(size(frameconditions));
for i = 1:numel(frameconditions),
  conditions{i} = '';
  if ~strcmpi(frameconditions{i},'any'),
    conditions{i} = frameconditions{i};
  end
  if ~strcmpi(flyconditions{i},'any'),
    if isempty(conditions{i}),
      conditions{i} = [conditions{i},'__'];
    end
    conditions{i} = [conditions{i},frameconditions{i}];
  end
end
  