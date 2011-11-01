% DebugSAGEExpDirConsistencyCheck

%% set up path

[~,computername] = system('hostname');
computername = strtrim(computername);

switch computername,
  
  case 'bransonk-lw1',
    
    addpath E:\Code\JCtrax\misc;
    addpath E:\Code\JCtrax\filehandling;
    addpath E:\Code\FlyBowlDataCapture;
    addpath('E:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
  
  case 'bransonk-lw2',

    addpath C:\Code\JCtrax\misc;
    addpath C:\Code\JCtrax\filehandling;
    addpath C:\Code\FlyBowlDataCapture;
    addpath('C:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    
  case 'bransonk-desktop',
    
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
    addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
    addpath /groups/branson/bransonlab/projects/olympiad/FlyBowlDataCapture;
    settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';

  otherwise
    
    warning('Unknown computer %s. Paths may not be setup correctly',computername);
    addpath C:\Code\JCtrax\misc;
    addpath C:\Code\JCtrax\filehandling;
    addpath C:\Code\FlyBowlDataCapture;
    addpath('C:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    
end

analysis_protocol = '20111005';
outfilename = 'SAGEExpDirConsistencyCheckResults20111021.txt';

global DISABLESAGEWAITBAR;
DISABLESAGEWAITBAR = true;

%% get all experiment directories

[expdirs,~,~,experiments] = getExperimentDirs('rootdir',rootdatadir,...
  'settingsdir',settingsdir,'protocol',analysis_protocol,...
  'daterange',{'20111008T000000','20111015T000000'});

[~,order] = sort({experiments.exp_datetime});

expdirs = expdirs(order(end:-1:1));
experiments = experiments(order(end:-1:1));

%% check each experiment

fid = fopen(outfilename,'w');
if fid < 0,
  error('Could not open file %s for output\n',outfilename);
end
isconsistent = nan(1,numel(expdirs));
filesmissing = nan(1,numel(expdirs));
for i = 1:numel(expdirs),
  expdir = experiments(i).file_system_path;
  fprintf('%s     ...    ',expdirs{i});
  [isconsistent(i),filesmissing(i)] = SAGEExpDirConsistencyCheck(expdir,...
    'analysis_protocol',analysis_protocol,...
    'settingsdir',settingsdir,...
    'docheckhist',false,...
    'outfilename',outfilename);
  if isconsistent(i),
    fprintf('consistent ');
  else
    fprintf('INCONSISTENT ');
  end
  if filesmissing(i),
    fprintf('FILESMISSING\n');
  else
    fprintf('filesexist\n');
  end
end

%% grab out only the inconsistent experiments

idx = find(~isconsistent);
fid = fopen(outfilename,'r');
fid2 = fopen('SAGEInconsisteExps20110907.txt','w');
for i = idx(:)',
  expdir = experiments(i).file_system_path;
  isfirst = true;
  didfind = false;
  while true,
    s = fgetl(fid);
    if ~ischar(s),
      if isfirst,
        fseek(fid,0,'bof');
        isfirst = false;
      else
        break;
      end
    end
    if strcmp(s,expdir),
      didfind = true;
      break;
    end
  end
  fprintf(fid2,'\n%s\n',expdir);
  while true,
    s = fgetl(fid);
    if ~ischar(s) || isempty(s),
      break;
    end
    fprintf(fid2,'%s\n',s);
  end
  
end
fclose(fid2);
fclose(fid);