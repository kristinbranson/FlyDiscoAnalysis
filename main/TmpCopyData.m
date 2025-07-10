function TmpCopyData(rootdir,outdir,rootdir0)

if nargin < 3,
  rootdir0 = rootdir;
end

l = mydir(rootdir,'isdir',true);
if isempty(l),
  return;
end

for i = 1:numel(l),
  
  l1 = mydir(fullfile(l{i},'score*mat'));
  if isempty(l1),
    TmpCopyData(l{i},fullfile(outdir,l{i}(numel(rootdir0)+1:end)));
    continue;
  end
  
  if ~exist(outdir,'dir'),
    mkdir(outdir);
  end
  for j = 1:numel(l1),
    fprintf('%s\n',l1{j});
    [~,name] = myfileparts(l1{j});
    unix(sprintf('cp %s %s',l1{j},fullfile(outdir,name)));
  end
  unix(sprintf('cp %s %s',fullfile(l{i},'trx.mat'),fullfile(outdir,'trx.mat')));
end