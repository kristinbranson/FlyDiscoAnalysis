function histperexp = MergeHistPerExp(histperexp,histperexpnew,expdir)

fns = fieldnames(histperexpnew);
for i = 1:numel(fns),
  fn = fns{i};
  nflies = numel(histperexpnew.(fn).Z);
  histperexpnew.(fn).expdir = repmat({expdir},[1,nflies]);
end


if isempty(histperexp),
  histperexp = histperexpnew;
  return;
end

for i = 1:numel(fns),
  
  fn = fns{i};
  
  % if this field is not yet set for histperexp, then just copy over
  if ~isfield(histperexp,fn),
    histperexp.(fn) = histperexpnew.(fn);
    continue;
  end
  
  % append
  fns2 = fieldnames(histperexpnew.(fn));
  for j = 1:numel(fns2),
    fn2 = fns2{j};
    histperexp.(fn).(fn2) = cat(1,histperexp.(fn).(fn2),histperexpnew.(fn).(fn2));
  end
  
end