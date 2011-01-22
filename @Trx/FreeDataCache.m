function FreeDataCache(obj,ndataadd)

while true,
  
  if obj.ndatacached + ndataadd <= obj.maxdatacached || obj.ndatacached == 0,
    break;
  end

  [~,j] = max(obj.perframehistory(:,2));
  fn = obj.perframehistory{j,1};
  for i = 1:obj.nexpdirs,
    if isfield(obj.datacached{i},fn),
      ndatacurr = 0;
      for j = 1:numel(obj.datacached{i}),
        ndatacurr = ndatacurr + numel(obj.datacached{i}(j).(fn));
      end
      obj.datacached{i} = rmfield(obj.datacached{i},fn);
    end
    obj.ndatacachedperexp(i) = obj.ndatacachedperexp(i) - ndatacurr;
    obj.ndatacached = obj.ndatacached - ndatacurr;
  end
  
  obj.perframehistory(j,:) = [];

end