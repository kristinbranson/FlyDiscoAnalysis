% change soft-linked perframe directories to real directories with
% soft-links to the mat files. 

rootoutputdir = '/groups/branson/bransonlab/projects/olympiad/HackHitData';

tmp = dir(rootoutputdir);
experiment_names = {tmp.name};
idx = cellfun(@isempty,regexp(experiment_names,'_\d{8}T\d{6}','once'));
experiment_names(idx) = [];
expdirs = cellfun(@(s) fullfile(rootoutputdir,s),experiment_names,'UniformOutput',false);

for i = 1:numel(expdirs),
  expdir = expdirs{i};
  perframedir = fullfile(expdir,'perframe');
  if ~exist(perframedir,'dir'),
    continue;
  end
  [~,link] = unix(sprintf('readlink %s',perframedir));
  link = strtrim(link);
  if isempty(link),
    continue;
  end
  fprintf('Fixing %s...\n',experiment_names{i});
  cmd = sprintf('rm %s',perframedir);
  unix(cmd);
  cmd = sprintf('mkdir %s',perframedir);
  unix(cmd);
  if ~exist(perframedir,'dir'),
    error('Could not create perframe directory %s',perframedir);
  end
  files = dir(link);
  file_names = setdiff({files.name},{'.','..'});
  for j = 1:numel(file_names),
    cmd = sprintf('ln -s %s %s',fullfile(link,file_names{j}),fullfile(perframedir,file_names{j}));
    unix(cmd);
    if ~exist(fullfile(perframedir,file_names{j}),'file'),
      error('Error creating link %s',fullfile(perframedir,file_names{j}));
    end
  end
end