function SymbolicCopyExampleExperiments(srcdir,destdir,debug)

% srcdir = '/groups/branson/bransonlab/taylora/flydisco/example-experiments';
% destdir = '/groups/branson/home/bransonk/behavioranalysis/code/example-experiments';
% 
% SymbolicCopyExampleExperiments(srcdir,destdir);

if nargin < 3,
  debug = false;
end

fs = dir(srcdir);
for i = 1:numel(fs),
  f = fs(i);
  if ismember(f.name,{'.','..'}),
    continue;
  end
  if f.isdir,
    srcdir1 = fullfile(srcdir,f.name);
    isexpdir = is_experiment_folder_given_contents(srcdir1);
    %[~,isexpdir] = parseExpDir(f.name);
    destdir1 = fullfile(destdir,f.name);
    if isexpdir,
      if exist(destdir1,'dir'),
        continue;
      end
      if debug,
        fprintf('SymbolicCopyExperimentDirectory(%s,%s)\n',srcdir1,destdir1);
      else
        SymbolicCopyExperimentDirectory(srcdir1,destdir1);
      end
    else
      if ~exist(destdir1,'dir'),
        if debug,
          fprintf('mkdir(%s)\n',destdir1);
        else
          mkdir(destdir1);
        end
      end
      SymbolicCopyExampleExperiments(srcdir1,destdir1,debug);
    end
  else
    assert(~f.isdir);
    if exist(fullfile(destdir,f.name),'file'),
      continue;
    end
    cmd = sprintf('ln -s %s %s',fullfile(srcdir,f.name),fullfile(destdir,f.name));
    if debug,
      fprintf('%s\n',cmd);
    else
      unix(cmd);
    end
  end
end