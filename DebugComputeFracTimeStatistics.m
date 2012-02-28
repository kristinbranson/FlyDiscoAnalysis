%% set up paths

[~,computername] = system('hostname');
computername = strtrim(computername);

switch computername,
  
  case 'bransonk-lw1',
    
    addpath E:\Code\JCtrax\misc;
    addpath E:\Code\JCtrax\filehandling;
    addpath('E:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
  
  case 'bransonk-lw2',

    addpath C:\Code\JCtrax\misc;
    addpath C:\Code\JCtrax\filehandling;
    addpath('C:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    
  case 'bransonk-desktop',
    
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
    addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
    settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';

  otherwise
    
    warning('Unknown computer %s. Paths may not be setup correctly',computername);
    addpath C:\Code\JCtrax\misc;
    addpath C:\Code\JCtrax\filehandling;
    addpath('C:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    
end

%% parameters

analysis_protocol = 'current';

tmpinputdir = '/groups/branson/home/kabram/clusterScripts/matlaboutWalk';

params = {'analysis_protocol',analysis_protocol,'settingsdir',settingsdir};

% SAGEParams = {'screen_type','non_olympiad_mushroombody','checkflags',true,'removemissingdata',true,...
%   'daterange',{'20110101T000000'},'rootdir',rootdatadir};
% SAGEParams = {'screen_type','primary','checkflags',true,'removemissingdata',true,...
%   'rootdir',rootdatadir};
SAGEParams = {'screen_type','primary','daterange',{'20120101T000000','20120129T000000'},...
  'checkflags',true,'removemissingdata',true};

ds = '20120203';

rootoutputdir = '/groups/branson/bransonlab/projects/olympiad/HackHitData';


%% get a list of experiment directories to process

data = SAGEListBowlExperiments(SAGEParams{:});

[~,order] = sort({data.exp_datetime});
order = order(end:-1:1);
data = data(order);

% filestoremove = {'*_labels.mat','scores_*.mat'};
% for i = 1:numel(data),
%   for j = 1:numel(filestoremove),
%     cmd = sprintf('rm %s',fullfile(data(i).file_system_path,filestoremove{j}));
%     unix(cmd);
%   end
% end
  

%% create output directories
for i = 1:numel(data),
  if mod(i,10) == 0, fprintf('%d / %d\n',i,numel(data)); end
  [~,experiment_name] = fileparts(data(i).file_system_path);
  inexpdir = fullfile(rootdatadir,experiment_name);
  outexpdir = fullfile(rootoutputdir,experiment_name);
  if exist(fullfile(outexpdir,'movie.ufmf'),'file'),
    continue;
  end
  SymbolicCopyExperimentDirectory(inexpdir,rootoutputdir);
  data(i).file_system_path = fullfile(rootoutputdir,experiment_name);
end


%% copy data from Mayank's directory

if ~isempty(tmpinputdir) && exist(tmpinputdir,'dir'),
  for i = 1:numel(data),
    if mod(i,100) == 0, fprintf('%d / %d\n',i,numel(data)); end
    [~,experiment_name] = fileparts(data(i).file_system_path);
    tmpexpdir = fullfile(tmpinputdir,experiment_name);
    if ~exist(tmpexpdir,'dir'),
      warning('Input directory %s does not exist',tmpexpdir);
      continue;
    end
    outexpdir = fullfile(rootoutputdir,experiment_name);
    if ~exist(outexpdir,'dir'),
      warning('Output directory %s does not exist',outexpdir);
      continue;
    end
    
    cmd = sprintf('cp %s/*.mat %s/.',tmpexpdir,outexpdir);
    unix(cmd);
  end
end

%% run in sequence

parfor i = 1:numel(data),
  [~,experiment_name] = fileparts(data(i).file_system_path);
  expdir = fullfile(rootoutputdir,experiment_name);
  [statsperfly,statsperexp] = ComputeFracTimeStatistics(expdir,params{:});
end

%% output to file to run in parallel

filename = sprintf('expdirs_fractime_%s.txt',ds);
fid = fopen(filename,'w');

for i = 1:numel(data),
  [~,experiment_name] = fileparts(data(i).file_system_path);
  tmpexpdir = fullfile(rootoutputdir,experiment_name);
  if exist(data(i).file_system_path,'dir') && ...
      exist(tmpexpdir,'dir'),
    fprintf(fid,'%s\n',tmpexpdir);
  end
end

fclose(fid);

%% see what failed!

filename = sprintf('expdirs_fractime_%s_rerun.txt',ds);
fid = fopen(filename,'w');
mindn = datenum(ds,'yyyymmdd');
behavior_names = {'stop','walk','chase','jump'};

for i = 1:numel(data),
  if mod(i,100) == 0, fprintf('%d / %d\n',i,numel(data)); end
  [~,experiment_name] = fileparts(data(i).file_system_path);
  tmpexpdir = fullfile(tmpinputdir,experiment_name);
  outexpdir = fullfile(rootoutputdir,experiment_name);
  if ~exist(outexpdir,'dir') || ~exist(tmpexpdir,'dir'),
    continue;
  end
  tmp = dir(fullfile(outexpdir,'fractime_stats.mat'));
  if isempty(tmp) || tmp.datenum < mindn,
    nbehaviorlabels = 0;
    for j = 1:numel(behavior_names),
      if exist(fullfile(tmpexpdir,[behavior_names{j},'_labels.mat']),'file'),
        nbehaviorlabels = nbehaviorlabels + 1;
      end
    end
    if nbehaviorlabels < numel(behavior_names),
      fprintf('Not enough labels for %s\n',experiment_name);
      continue;
    end
    fprintf(fid,'%s\n',outexpdir);
  end
end
fclose(fid);

%% 

filename = 'expdirs_hackhits_20111215_archived.txt';
fid = fopen(filename,'w');
for i = 1:numel(probarchived),
  fprintf(fid,'%s\n',fullfile(rootdatadir,probarchived{i}));
end
fclose(fid);