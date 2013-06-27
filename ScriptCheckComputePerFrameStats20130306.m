mindatenum = datenum('20130317','yyyymmdd');
isok = false(1,numel(expdirs));
for i = find(~isok),
  expdir = expdirs{i};
  statsfile = fullfile(expdir,'stats_perframe.mat');
  tmp = dir(statsfile);
  if isempty(tmp),
    fprintf('%s does not exist\n',statsfile);
    continue;
  elseif tmp(1).datenum < mindatenum,
    fprintf('%s is too old (%s)\n',statsfile,tmp(1).date);
    continue;
  end
  
  histfile = fullfile(expdir,'hist_perframe.mat');
  tmp = dir(histfile);
  if isempty(tmp),
    fprintf('%s does not exist\n',histfile);
    continue;
  elseif tmp(1).datenum < mindatenum,
    fprintf('%s is too old (%s)\n',histfile,tmp(1).date);
    continue;
  end
  
  
  isok(i) = true;
  
end