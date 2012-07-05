function SymbolicCopyFlyBowlDataCaptureFiles(expdir,rootoutputdir,varargin)

filestocopy = {'FlyBowlDataCaptureParams_EP*',...
  'Log.txt','Metadata*','movie.ufmf','QuickStats.png','Quickstats.txt',...
  'SUCCESS','temperature.txt','ufmf_diagnostics.txt','ufmf_log.txt'};

[filestocopy,dosoftlink] = myparse(varargin,'filestocopy',filestocopy,'dosoftlink',true);

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

doremove = false(1,numel(filestocopy));
for i = 1:numel(filestocopy),
  infile = fullfile(expdir,filestocopy{i});
  if any(filestocopy{i} == '*'),
    files = dir(infile);
    if isempty(files),
      doremove(i) = true;
    else
      filestocopy{i} = files(1).name;
      for j = 2:numel(files),
        filestocopy{end+1} = files(j).name; %#ok<AGROW>
      end
    end
  end
end
filestocopy(doremove) = [];

for i = 1:numel(filestocopy),
  infile = fullfile(expdir,filestocopy{i});
  if ~exist(infile,'file'),
    continue;
  end
  outfile = fullfile(outexpdir,filestocopy{i});
  if ~exist(outfile,'file'),
    if dosoftlink,
      cmd = sprintf('ln -s %s %s',infile,outfile);
    else
      cmd = sprintf('cp %s %s',infile,outfile);
    end
    fprintf('%s -> %s\n',infile,outfile);
    system(cmd);
    if ~exist(outfile,'file'),
      error('File %s not created successfully',outfile);
    end
  end
end
