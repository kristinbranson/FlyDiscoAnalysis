function [colors,markers] = SelectPlotParams(plotparams,frameconditions,flyconditions)

nfns = numel(frameconditions);
colors = repmat({[0,0,0]},[1,nfns]);
markers = repmat({'o'},[1,nfns]);
for fni = 1:nfns,
  framematch = strcmp(plotparams.frameconditions,frameconditions{fni});
  flymatch = strcmp(plotparams.flyconditions,flyconditions{fni});
  i = find(framematch & flymatch,1);
  if isempty(i),
    i = find(framematch&strcmp(plotparams.flyconditions,'*'),1);
    if isempty(i),
      i = find(flymatch&strcmp(plotparams.frameconditions,'*'),1);
    end
    if isempty(i),
      continue;
    end
  end
  colors{fni} = plotparams.colors{i};
  markers{fni} = plotparams.markers{i};  
end