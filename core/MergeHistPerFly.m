function histperfly = MergeHistPerFly(histperfly,histperflynew,expdir)

fns = fieldnames(histperflynew);
for i = 1:numel(fns),
  fn = fns{i};
  nflies = numel(histperflynew.(fn).Z);
  histperflynew.(fn).expdir = repmat({expdir},[1,nflies]);
end


if isempty(histperfly),
  histperfly = histperflynew;
  return;
end

for i = 1:numel(fns),
  
  fn = fns{i};
  
  % if this field is not yet set for histperfly, then just copy over
  if ~isfield(histperfly,fn),
    histperfly.(fn) = histperflynew.(fn);
    continue;
  end
  
  % append
  fns2 = fieldnames(histperflynew.(fn));
  for j = 1:numel(fns2),
    fn2 = fns2{j};
    histperfly.(fn).(fn2) = cat(2,histperfly.(fn).(fn2),histperflynew.(fn).(fn2));
  end
  
end