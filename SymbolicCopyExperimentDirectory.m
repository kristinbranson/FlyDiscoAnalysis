function SymbolicCopyExperimentDirectory(expdir,rootoutputdir)

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
  infile = fullfile(expdir,subfiles{i});
  outfile = fullfile(outexpdir,subfiles{i});
  if ~exist(outfile,'file'),
    cmd = sprintf('ln -s %s %s',infile,outfile);
    system(cmd);
    if ~exist(outfile,'file'),
      error('Soft-linked file %s not created successfully',outfile);
    end
  end
end

inperframedir = fullfile(expdir,'perframe');
outperframedir = fullfile(outexpdir,'perframe');
if exist(inperframedir,'dir'),
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
      cmd = sprintf('ln -s %s %s',infile,outfile);
      system(cmd);
      if ~exist(outfile,'file'),
        error('Soft-linked file %s not created successfully',outfile);
      end
    end
  end
end