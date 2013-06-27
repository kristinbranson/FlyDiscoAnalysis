function [expidxsample,n] = CreateControlDataSample(nexpsperset_line,nexpsperset_control,set2exps_control)

nexps_line = sum(nexpsperset_line);
nsets_control = numel(nexpsperset_control);
maxnexps = max(nexpsperset_line);
counts = hist(nexpsperset_line,1:maxnexps);

expidxsample = nan(1,nexps_line);
n = nan(1,nexps_line);
isleft = true(1,nsets_control);

off = 0;
for minnexps = maxnexps:-1:1,
  
  % choose some sets
  isallowed = nexpsperset_control >= minnexps;
  idxsample = randsample(find(isleft&isallowed),counts(minnexps),false);
  isleft(idxsample) = false;
  
  % choose some experiments
  n(off+1:end+counts(minnexps)*minnexps) = minnexps;
  for seti = 1:counts(minnexps),
    expidxsample(off+1:off+minnexps) = set2exps_control{idxsample(seti)}(randsample(nexpsperset_control(idxsample(seti)),minnexps));
    off = off + minnexps;
  end
  
end