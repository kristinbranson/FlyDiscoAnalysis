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
rootdatadir = '/groups/branson/bransonlab/projects/olympiad/HackHitData';

params = {'analysis_protocol',analysis_protocol};
%  'rootoutputdir',rootoutputdir};

% SAGEParams = {'screen_type','primary','checkflags',true,'removemissingdata',true,...
%   'daterange',{'20110301T000000'}};
SAGEParams = {'screen_type','non_olympiad_mushroombody','checkflags',true,'removemissingdata',true,...
  'daterange',{'20110101T000000'},'rootdir',rootdatadir};


%% get a list of experiment directories to process

data = SAGEListBowlExperiments(SAGEParams{:});

[~,order] = sort({data.exp_datetime});
order = order(end:-1:1);
data = data(order);

%% run in sequence

for i = 1:numel(data),
  fprintf('%d / %d\n',i,numel(data));
  expdir = data(i).file_system_path;
  [statsperfly,statsperexp] = ComputeHackHitCategoryStatistics_KB(expdir,params{:});
end


%%

nperfile = 1000;

for filei = 1:ceil(numel(data)/nperfile),

  filename = sprintf('expdirs_hackhits_20111215_%02d.txt',filei);
  fid = fopen(filename,'w');

  for i = nperfile*(filei-1)+1:min(numel(data),nperfile*filei),
    fprintf(fid,'%s\n',data(i).file_system_path);
  end

  fclose(fid);
end

%%

filename = 'expdirs_hackhits_20111215_02.txt';
fid = fopen(filename,'w');
off = 1000;

for i = off+1:numel(data),
  fprintf(fid,'%s\n',data(i).file_system_path);
end

fclose(fid);

%% 

tmp = dir(rootoutputdir);
probarchived = {};
for i = 1:numel(tmp),
  if ismember(tmp(i).name,{'.','..'}),
    continue;
  end
  fn = fullfile(rootoutputdir,tmp(i).name,'hackhitstats.mat');
  if ~exist(fn,'file'),
    fprintf('Missing hackhitstats.mat for %s\n',tmp(i).name);
    probarchived{end+1} = tmp(i).name;
  end
end

%% 

filename = 'expdirs_hackhits_20111215_archived.txt';
fid = fopen(filename,'w');
for i = 1:numel(probarchived),
  fprintf(fid,'%s\n',fullfile(rootdatadir,probarchived{i}));
end
fclose(fid);