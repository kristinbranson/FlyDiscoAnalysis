function outexpdir = SymbolicCopyExperimentDirectory(expdir,rootoutputdir,varargin)

[copyfiles,ignorefiles,dosoftlink] = myparse(varargin,...
  'copyfiles',{},'ignorefiles',{},'dosoftlink',true);

% create the root directory
if ~exist(rootoutputdir,'dir'),
  [success,msg] = mkdir(rootoutputdir);
  if ~success,
    error('Could not create root output directory %s: %s',rootoutputdir,msg);
  end
end
  
% create the experiment directory
[~,basename] = fileparts(expdir);
outexpdir = fullfile(rootoutputdir,basename);
if ~exist(outexpdir,'dir'),
  [success,msg] = mkdir(outexpdir);
  if ~success,
    error('Could not create output directory %s: %s',outexpdir,msg);
  end
end

subfiles = dir(expdir);
subfiles = setdiff({subfiles.name},{'.','..','perframe'});

for i = 1:numel(subfiles),
  
  if ~isempty(copyfiles) && ~ismember(subfiles{i},copyfiles),
    continue;
  end
  
  if ismember(subfiles{i},ignorefiles),
    continue;
  end
  
  infile = fullfile(expdir,subfiles{i});
  outfile = fullfile(outexpdir,subfiles{i});
  if ~exist(outfile,'file'),
    if dosoftlink,
      [~,link] = unix(sprintf('readlink %s',infile));
      link = strtrim(link);
      if ~isempty(link),        
        cmd = sprintf('ln -s %s %s',link,outfile);
      else
        cmd = sprintf('ln -s %s %s',infile,outfile);
      end
    else
      cmd = sprintf('cp -r %s %s',infile,outfile);
    end
    system(cmd);
    if ~exist(outfile,'file'),
      error('Soft-linked file %s not created successfully',outfile);
    end
  end
end

inperframedir = fullfile(expdir,'perframe');
outperframedir = fullfile(outexpdir,'perframe');
if exist(inperframedir,'dir') && (isempty(copyfiles) || ismember('perframe',copyfiles)) && ~ismember('perframe',ignorefiles),
  if ~exist(outperframedir,'dir'),
    cmd = sprintf('mkdir %s',outperframedir);
    unix(cmd);
  end
  subfiles = dir(inperframedir);
  subfiles = setdiff({subfiles.name},{'.','..','a_mm.mat','b_mm.mat','x_mm.mat','y_mm.mat','theta_mm.mat'});

  for i = 1:numel(subfiles),
    infile = fullfile(inperframedir,subfiles{i});
    outfile = fullfile(outperframedir,subfiles{i});
    if ~exist(outfile,'file'),
      if dosoftlink,
        cmd = sprintf('ln -s %s %s',infile,outfile);
      else
        cmd = sprintf('cp %s %s',infile,outfile);
      end
      system(cmd);
      if ~exist(outfile,'file'),
        error('Soft-linked file %s not created successfully',outfile);
      end
    end
  end
end