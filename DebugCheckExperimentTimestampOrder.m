%% set up path

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;

settinsgdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
datalocparamsfilestr = 'dataloc_params.txt';

%% parameters

analysis_protocol = 'current';

%% list all experiments

allexps = SAGEListBowlExperiments('flag_aborted','0','checkflags',false,'screen_type','primary');
expdirs = {allexps.file_system_path};

%% allocate

filesmissing = cell(1,numel(allexps));
datalocparams = ReadParams(fullfile(settingsdir,analysis_protocol,datalocparamsfilestr));
perframefnsfile = fullfile(settingsdir,analysis_protocol,datalocparams.perframefnsfilestr);
allperframefns = importdata(perframefnsfile);

%% run once to get rules

i = 1;
[isbroken1,filesmissing{i},filetypes,rules] = CheckExperimentTimestampOrder(expdirs{i},'debug',false,'outputdir','','allperframefns',allperframefns);
nrules = size(rules,1);

isbroken = nan(numel(allexps),nrules);
isbroken(i,:) = isbroken1;

%% run lots of times

%order = randperm(numel(allexps));

parfor i = 1:numel(allexps),
  %i = order(ii);
  expdir = expdirs{i};
  [isbroken(i,:),filesmissing{i}] = CheckExperimentTimestampOrder(expdirs{i},'debug',false,'outputdir','','allperframefns',allperframefns);
  if mod(i,100) == 0,
    fprintf('Processed %d experiments / %d\n',i,numel(allexps));
%     fprintf('Frac is broken: ');
%     fprintf('%f ',sum(isbroken(order(1:ii),:),1)/i);
%     fprintf('\n');
  end
end