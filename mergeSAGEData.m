function data1 = mergeSAGEData(data1,data2)

expnames1 = {data1.experiment_name};
expnames2 = {data2.experiment_name};
fns1 = fieldnames(data1);
fns2 = fieldnames(data2);
[isold,idx] = ismember(expnames2,expnames1);
idxold = find(isold);
idxnew = find(~isold);
idxnotnew = find(~ismember(expnames1,expnames2));
fns2not1 = setdiff(fns2,fns1);
fns1not2 = setdiff(fns1,fns2);
for i2 = idxold(:)',
  i1 = idx(i2);
  for j = 1:numel(fns2not1),
    fn = fns2not1{j};
    data1(i1).(fn) = data2(i2).(fn);
  end
end
for i2 = idxnew(:)',
  i1 = numel(data1)+1;
  for j = 1:numel(fns2),
    fn = fns2{j};
    data1(i1).(fn) = data2(i2).(fn);
  end
end
if ~isempty(idxnew) && ~isempty(fns1not2),
  fprintf('There are experiments in data2 not in data1, and there will be missing data for the fields:\n');
  fprintf('%s\n',fns1not2{:});
  fprintf('Experiments that are new:\n');
  fprintf('%s\n',expnames2{idxnew});
end
if ~isempty(idxnotnew) && ~isempty(fns2not1),
  fprintf('There are experiments in data1 not in data2, and there will be missing data for the fields:\n');
  fprintf('%s\n',fns2not1{:});
  fprintf('Experiments that are not in 2:\n');
  fprintf('%s\n',expnames1{idxnotnew});
end
  