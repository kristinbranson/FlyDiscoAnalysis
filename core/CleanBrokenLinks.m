function CleanBrokenLinks(dirname)

files = dir(dirname);
for i = 1:numel(files),
  filename = fullfile(dirname,files(i).name);
  if ~exist(filename,'file'),
    fprintf('Deleting %s\n',filename);
    unix(sprintf('rm %s',filename));
  end
end