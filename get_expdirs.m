function [expdirs,expdir_parts] = get_expdirs(rootdir)

allfiles = dir(rootdir);
allfiles = allfiles([allfiles.isdir]);
expdirs = {};
expdir_parts = [];
for i = 1:length(allfiles),
  filecurr = allfiles(i).name;
  [partscurr,isexpdir] = parseExpDir(filecurr);
  if isexpdir,
    expdirs{end+1} = filecurr; %#ok<*AGROW>
    if isempty(expdir_parts),
      expdir_parts = partscurr;
    else
      expdir_parts(end+1) = partscurr;
    end
  end
end