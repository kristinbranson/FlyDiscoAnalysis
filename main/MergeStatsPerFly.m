function statsperfly = MergeStatsPerFly(statsperfly,statsperflynew,expdir)

fns = fieldnames(statsperflynew);
for i = 1:numel(fns),
  fn = fns{i};
  nflies = numel(statsperflynew.(fn).Z);
  statsperflynew.(fn).expdir = repmat({expdir},[1,nflies]);
end


if isempty(statsperfly),
  statsperfly = statsperflynew;
  return;
end

for i = 1:numel(fns),
  
  fn = fns{i};
  
  % if this field is not yet set for statsperfly, then just copy over
  if ~isfield(statsperfly,fn),
    statsperfly.(fn) = statsperflynew.(fn);
    continue;
  end
  
  % append
  fns2 = fieldnames(statsperflynew.(fn));
  for j = 1:numel(fns2),
    fn2 = fns2{j};
    statsperfly.(fn).(fn2) = cat(2,statsperfly.(fn).(fn2),statsperflynew.(fn).(fn2));
  end
  
end