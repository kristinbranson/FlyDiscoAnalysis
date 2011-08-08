% check directory-SAGE consistency

%% set up paths
[~,computername] = system('hostname');
computername = strtrim(computername);

switch computername,
  
  case 'bransonk-lw1',

    addpath E:\Code\SAGE\MATLABInterface\Trunk;
    addpath E:\Code\JCtrax\misc;
    rootdir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
    
  case 'bransonk-lw2',

    addpath C:\Code\SAGE\MATLABInterface\Trunk;
    addpath C:\Code\JCtrax\misc;
    rootdir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';

  case 'bransonk-desktop',
    
    addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
    rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
    settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    
  case 'robiea-ww1'
    
    addpath C:\Users\robiea\Documents\Code_versioned\SAGE\MATLABInterface\Trunk;
    addpath C:\Users\robiea\Documents\Code_versioned\JCtrax\misc;
    rootdir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    settingsdir = 'C:\Users\robiea\Documents\Code_versioned\FlyBowlAnalysis\settings';
    
  case 'robiea-ws'
    addpath '/groups/branson/home/robiea/Code_versioned/SAGE/MATLABInterface/Trunk';
    addpath '/groups/branson/home/robiea/Code_versioned/JCtrax/misc';
    rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
    settingsdir = '/groups/branson/home/robiea/Code_versioned/FlyBowlAnalysis/settings';
    
  otherwise
    
    warning('Unknown computer %s. Paths may not be setup correctly',computername);
    addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
    rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
    settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    
end

%% get all experiment directories

dirs = dir(fullfile(rootdir,'*_*'));
dirnames = {dirs([dirs.isdir]).name};
issuccess = false(1,numel(dirnames));
for i = 1:numel(dirnames),
  issuccess(i) = exist(fullfile(rootdir,dirnames{i},'movie.ufmf'),'file') && exist(fullfile(rootdir,dirnames{i},'Metadata.xml'),'file');
end
fprintf('%d experiment directories found, %d contain file Metadata.xml and movie.ufmf\n',numel(dirnames),nnz(issuccess));
dirnames = dirnames(issuccess);

%% get all experiments in SAGE

allexps = SAGEListBowlExperiments('checkflags',false,'removemissingdata',false,'rootdir',rootdir);
experiment_names = {allexps.experiment_name};

%% compare

dir_experiment_names = cellfun(@(s) ['FlyBowl_',s],dirnames,'UniformOutput',false);
dir_but_not_SAGE = setdiff(dir_experiment_names,experiment_names);
SAGE_but_not_dir = setdiff(experiment_names,dir_experiment_names);

%% sort

datetime = cell(1,numel(dir_but_not_SAGE));
for i = 1:numel(dir_but_not_SAGE),
  res = parseExpDir(dir_but_not_SAGE{i});
  datetime{i} = res.date;
end
  
[~,order] = sort(datetime);
dir_but_not_SAGE = dir_but_not_SAGE(order);

datetime = cell(1,numel(SAGE_but_not_dir));
for i = 1:numel(SAGE_but_not_dir),
  res = parseExpDir(SAGE_but_not_dir{i});
  datetime{i} = res.date;
end
  
[~,order] = sort(datetime);
SAGE_but_not_dir = SAGE_but_not_dir(order);

%% look for aborted

isaborted = false(1,numel(dir_but_not_SAGE));
for i = 1:numel(dir_but_not_SAGE),
  isaborted(i) = exist(fullfile(rootdir,dir_but_not_SAGE{i}(9:end),'ABORTED'),'file');
end



%% explain

reasons = cell(1,numel(dir_but_not_SAGE));
for i = 1:numel(dir_but_not_SAGE),
  if isaborted(i),
    reasons{i} = 'aborted';
    continue;
  end
  res = parseExpDir(dir_but_not_SAGE{i});
  if numel(res.date) >= 4 && strcmp(res.date(1:4),'2010'),
    reasons{i} = 'old';
    continue;
  end
  if isempty(dir(fullfile(rootdir,dir_but_not_SAGE{i}(9:end),'ctrax_results_movie*.avi'))),
    reasons{i} = 'analysis pipeline failed';
    continue;
  end
  
  %disp(dir_but_not_SAGE{i});
  %dir(fullfile(rootdir,dir_but_not_SAGE{i}(9:end)))
  %reasons{i} = input('\nReason: ','s');
  
end

%% print results

fid = fopen(sprintf('CheckFileSystemSAGEConsistency_%s.txt',datestr(now,'yyyymmddTHHMMSS')),'w');
fprintf(fid,'%d experiments found in directory structure but not in SAGE flattened view:\n',numel(dir_but_not_SAGE));
for i = 1:numel(dir_but_not_SAGE),
  fprintf(fid,'  %s,%s\n',dir_but_not_SAGE{i},reasons{i});
end
fprintf(fid,'\n\n\n%d experiments found in SAGE flattened view but not in directory structure:\n',numel(SAGE_but_not_dir));
fprintf(fid,'  %s\n',SAGE_but_not_dir{:});
fclose(fid);
