% check quarantine

rootdir = 'O:\Olympiad_Screen\fly_bowl';
required_files = {'movie.ufmf','Metadata.xml'};

%%

quarantine_dirs = dir(fullfile(rootdir,'*quarantine*'));
quarantine_dirs = {quarantine_dirs([quarantine_dirs.isdir]).name};

pipeline_dirs = dir(rootdir);
pipeline_dirs = {pipeline_dirs([pipeline_dirs.isdir]).name};
pipeline_dirs = pipeline_dirs(~cellfun(@isempty,regexp(pipeline_dirs,'^0\d','once')));
pipeline_dirs = setdiff(pipeline_dirs,quarantine_dirs);

for i = 1:numel(quarantine_dirs),
  
  expdirs = dir(fullfile(rootdir,quarantine_dirs{i},'*_*'));
  expdirs = {expdirs([expdirs.isdir]).name};
  fprintf('%s:\n',quarantine_dirs{i});
  for j = 1:numel(expdirs),
    expdir = fullfile(rootdir,quarantine_dirs{i},expdirs{j});
    issuccess = exist(fullfile(expdir,'SUCCESS'),'file');
    isaborted = exist(fullfile(expdir,'ABORTED'),'file');
    isnotstarted = ~isempty(regexp(expdirs{j},'notstarted','once'));
    if isaborted || isnotstarted || ~issuccess,
      continue;
    end
    if any(cellfun(@(s) ~exist(fullfile(expdir,s),'file'),required_files)),
      continue;
    end
    fprintf('  %s\n',expdirs{j});
  end
  
end

for i = 1:numel(pipeline_dirs)-1,
  
  expdirs = dir(fullfile(rootdir,pipeline_dirs{i},'*_*'));
  expdirs = {expdirs([expdirs.isdir]).name};
  fprintf('%s:\n',pipeline_dirs{i});
  for j = 1:numel(expdirs),
    expdir = fullfile(rootdir,pipeline_dirs{i},expdirs{j});
    issuccess = exist(fullfile(expdir,'SUCCESS'),'file');
    isaborted = exist(fullfile(expdir,'ABORTED'),'file');
    isnotstarted = ~isempty(regexp(expdirs{j},'notstarted','once'));
    if isaborted || isnotstarted || ~issuccess,
      continue;
    end
    if any(cellfun(@(s) ~exist(fullfile(expdir,s),'file'),required_files)),
      continue;
    end
    fprintf('  %s\n',expdirs{j});
  end
  
end
