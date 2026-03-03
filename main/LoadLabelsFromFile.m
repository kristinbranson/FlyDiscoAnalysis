function labelidx = LoadLabelsFromFile(labelfilename,nframes,firstframes)

nflies = numel(nframes);

if ~exist(labelfilename,'file'),
  error('File %s does not exist',labelfilename);
end

loadedlabels = load(labelfilename,'t0s','t1s','names','flies','off','timestamp');
labelidx = cell(1,nflies);
for i = 1:numel(loadedlabels.t0s),
  fly = loadedlabels.flies(:,i);
  labelidx{fly} = false(1,nframes(fly));
  off = loadedlabels.off(i);
  t0s = loadedlabels.t0s{i} - off;
  i0s = t0s + 1 - firstframes(fly);
  t1s = loadedlabels.t1s{i} - off;
  i1s = t1s + 1 - firstframes(fly);
  for j = 1:numel(i0s),
    labelidx{fly}(i0s(j):i1s(j)-1) = true;
  end
  
end
