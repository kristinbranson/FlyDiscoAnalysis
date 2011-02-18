% CombineLineDataScript

% set up path
if ispc,
  addpath E:\Code\JCtrax\misc;
  addpath E:\Code\JCtrax\filehandling;
else
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
  addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
end
 
%% locations of data

settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
analysis_protocol = '20110211';
datalocparamsfilestr = 'dataloc_params.txt';
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
timestamp = datestr(now,30);

%% get all linenames

[expdirs,expdir_reads,~,experiments,rootdir,~] = ...
  getExperimentDirs('protocol',analysis_protocol,'subreadfiles',...
  {dataloc_params.statsperframematfilestr,dataloc_params.histperframematfilestr});

linenames = unique({experiments.line});

%% loop through lines

for i = 1:numel(linenames),
  
  linename = linenames{i};
  
  fprintf('Combining data for line %s...\n',linename);
  
  % directory containing the data for this line
  linedatadircurr = fullfile(dataloc_params.linedatadir,linename);
  if ~exist(linedatadircurr,'file'),
    mkdir(linedatadircurr);
  end
  
  % directory for this instance of data
  outdir = fullfile(linedatadircurr,timestamp);
  if ~exist(outdir,'file'),
    mkdir(outdir);
  end

  FlyBowlCombineExperiments(rootdir,outdir,...
    'settingsdir',settingsdir,...
    'analysis_protocol',analysis_protocol,...
    'datalocparamsfilestr',datalocparamsfilestr,...
    'linename',linename,...
    'plottitle',linename);
  
end