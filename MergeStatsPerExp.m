function statsperexp = MergeStatsPerExp(statsperexp,statsperexpnew,expdir)

fns = fieldnames(statsperexpnew);
for i = 1:numel(fns),
  fn = fns{i};
  nflies = numel(statsperexpnew.(fn).Z);
  statsperexpnew.(fn).expdir = repmat({expdir},[1,nflies]);
end


if isempty(statsperexp),
  statsperexp = statsperexpnew;
  return;
end

for i = 1:numel(fns),
  
  fn = fns{i};
  
  % if this field is not yet set for statsperexp, then just copy over
  if ~isfield(statsperexp,fn),
    statsperexp.(fn) = statsperexpnew.(fn);
    continue;
  end
  
  % append
  fns2 = fieldnames(statsperexpnew.(fn));
  for j = 1:numel(fns2),
    fn2 = fns2{j};
    statsperexp.(fn).(fn2) = cat(1,statsperexp.(fn).(fn2),statsperexpnew.(fn).(fn2));
  end
  
end