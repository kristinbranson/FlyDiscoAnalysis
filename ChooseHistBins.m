%% set up path
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
realrootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
rootdir = '/groups/branson/bransonlab/projects/olympiad/HackHitData';
rootdir_fixed = '';

analysis_protocol = '20120330';
datalocparamsfilestr = 'dataloc_params.txt';
params = {'analysis_protocol',analysis_protocol,...
  'datalocparamsfilestr',datalocparamsfilestr,...
  'settingsdir',settingsdir};

behaviornames = {'chase','walk','stop','jump'};

fps = 30.3445;

%% parameters

datalocparamsfile = fullfile(settingsdir,analysis_protocol,'dataloc_params.txt');
dataloc_params = ReadParams(datalocparamsfile);

%% choose experiments to histogram

moviefilestr = 'movie.ufmf';
trxfilestr = dataloc_params.ctraxfilestr;
annfilestr = 'movie.ufmf.ann';
% 
% expfiles = dir(fullfile(rootoutputdir,'*_*'));
% experiments_all = [];
% for i = 1:numel(expfiles),
%   expdir = fullfile(rootoutputdir,expfiles(i).name);
%   [res,success] = parseExpDir(expdir);
%   if success,
%     res.file_system_path = expdir;
%     res.exp_datetime = res.date;
%     if isempty(experiments_all),
%       experiments_all = res;
%     else
%       experiments_all(end+1) = res; %#ok<SAGROW>
%     end
%   end
% end

experiment_params = struct;
% root data dir
experiment_params.rootdir = rootdir;
% what dates should we analyze
%experiment_params.daterange = {'20110201T000000','20111221T000000'};
% remove failures
experiment_params.checkflags = false;
% remove missing data
experiment_params.removemissingdata = false;
% azanchir experiments
experiment_params.screen_type = 'primary';
% how to weight various things
weight_order = {'genotype','date','rig','bowl'};
% required files
subreadfiles = {'chase_labels.mat','stop_labels.mat','walk_labels.mat','jump_labels.mat'};
% what lines
%experiment_params.linename = '';
% 
ngal4 = 100;
ncontrols = 100;
tmpparams = struct2paramscell(experiment_params);

%[~,~,~,experiments_all] = getExperimentDirs(tmpparams{:});
experiments_all = SAGEListBowlExperiments(tmpparams{:});

baddata = false(1,numel(experiments_all));
for i = 1:numel(experiments_all),
  if ~exist(experiments_all(i).file_system_path,'dir'),
    baddata(i) = true;
    continue;
  end
  [~,experiment_name] = myfileparts(experiments_all(i).file_system_path);
%   annfile = fullfile(rootdir,experiment_name,annfilestr);
%   if ~exist(annfile,'file'),
%     baddata(i) = true;
%   end
  for j = 1:numel(subreadfiles),
    if ~exist(fullfile(rootdir,experiment_name,subreadfiles{j}),'file'),
      baddata(i) = true;
    end
  end
end
experiments_all(baddata) = [];

%experiments_all = rmfield(experiments_all,'date');

[experiments,ngal4_chosen,ncontrols_chosen] = choose_expdirs(experiments_all,ngal4,ncontrols,'weight_order',weight_order);

expdirs = {experiments.file_system_path};
nexpdirs = numel(expdirs);



%% load data

% obj = Trx('settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
% for i = 1:numel(expdirs),
%   obj.AddExpDir(expdirs{i});
% end

%% read parameters

perframefnsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.perframefnsfilestr);
perframefns = importdata(perframefnsfile);
prctile_store = 10;

%histparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histparamsfilestr);

histperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframefeaturesfilestr);
hist_perframefeatures = ReadHistPerFrameFeatures2(histperframefeaturesfile);
hist_params = ReadParams(histparamsfile);

%% loop over per-frame parameters

% data structure we will store everything in
bins = struct;
perframefns = setdiff({hist_perframefeatures.field},{'duration'});

% add in behavior labels
for i = 1:numel(behaviornames),
  perframefns{end+1} = ['duration_',behaviornames{i}];
end
  
for fni = 1:numel(perframefns),
  
  fn = perframefns{fni};
  
  if strcmp(fn,'sex'),
    continue;
  end
  
  fprintf('Choosing bins for %s...\n',fn);
  bins.(fn) = struct;
  
  % store smallest and largest data to get lims for histogramming
  datasmall = [];
  datalarge = [];

  fprintf('Choosing lims...\n');
  % loop through experiments
  for expi = 1:numel(expdirs),
    
    expdir = expdirs{expi};
    %fprintf('exp(%d) = %s\n',expi,expdir);

    m = regexp(fn,'^duration_(.+)$','tokens','once');
    if ~isempty(m),
%       trxfilename = fullfile(expdir,dataloc_params.trxfilestr);
%       trxdata = load(trxfilename);
%       firstframes = [trxdata.trx.firstframe];
%       nframes = [trxdata.trx.nframes];

      filename = fullfile(expdir,[m{1},'_labels.mat']);
      labeldata = load(filename);
      datacurr = (cell2mat(labeldata.t1s)-cell2mat(labeldata.t0s))/fps;
    else
      filename = fullfile(expdir,dataloc_params.perframedir,[fn,'.mat']);
      if ~exist(filename,'file'), 
        warning('Could not find file %s, skipping',filename);
        continue;
      end
      % get per-frame data for these flies
      tmp = load(filename,'data');
      datacurr = cell2mat(tmp.data);
    end
    thresh = prctile(datacurr,[prctile_store,100-prctile_store]);
    datasmall = [datasmall,datacurr(datacurr <= thresh(1))]; %#ok<AGROW>
    datalarge = [datalarge,datacurr(datacurr >= thresh(2))]; %#ok<AGROW>
    
  end
  
  % compute actual lim percentiles from the stored boundary data
  bins.(fn).lims = [nan,nan];

  mindata = min(datasmall);
  maxdata = max(datalarge);
  % set lower bound to 0 if this makes sense
  if mindata >= 0 && mindata / maxdata < 1e-3,
    bins.(fn).lims(1) = 0;
    fprintf('Choosing lb = 0 for %s\n',fn);
  % set lower bound to -pi if this makes sense
  elseif maxdata > -pi && abs(mindata + pi) / (maxdata + pi) < 1e-3,
    bins.(fn).lims(1) = -pi;
    fprintf('Choosing lb = -pi for %s\n',fn);
  % set lower bound to -pi/2 if this makes sense
  elseif maxdata > -pi/2 && abs(mindata + pi/2) / (maxdata + pi/2) < 1e-3,
    bins.(fn).lims(1) = -pi/2;
    fprintf('Choosing lb = -pi/2 for %s\n',fn);
  % set lower bound to -1 if this makes sense
  elseif maxdata > -1 && abs(mindata + 1) / (maxdata + 1) < 1e-3,
    bins.(fn).lims(1) = -1;
    fprintf('Choosing lb = -1 for %s\n',fn);
  else
    bins.(fn).lims(1) = prctile(datasmall,hist_params.lim_prctile(1)*prctile_store);
    fprintf('Choosing lb from percentile for %s\n',fn);
  end
  
  % set upper bound to pi if this makes sense
  if mindata < pi && abs(maxdata - pi) / (pi - mindata) < 1e-3,
    bins.(fn).lims(2) = pi;
    fprintf('Choosing ub = pi for %s\n',fn);
  % set upper bound to pi/2 if this makes sense
  elseif mindata < pi/2 && abs(maxdata - pi/2) / (pi/2 - mindata) < 1e-3,
    bins.(fn).lims(2) = pi/2;
    fprintf('Choosing ub = pi/2 for %s\n',fn);
  % set upper bound to 1 if this makes sense
  elseif mindata < 1 && abs(maxdata - 1) / (1 - mindata) < 1e-3,
    bins.(fn).lims(2) = 1;
    fprintf('Choosing ub = 1 for %s\n',fn);
  else    
    bins.(fn).lims(2) = prctile(datalarge,100-(100-hist_params.lim_prctile(2))*prctile_store);
    fprintf('Choosing ub from percentile for %s\n',fn);
  end

  % edges based on linear spacing
  [bins.(fn).edges_linear,bins.(fn).centers_linear] = ...
    SelectHistEdges(hist_params.nbins,bins.(fn).lims,'linear');
  
  % edges based on log spacing
  if mindata < 0,
    [bins.(fn).edges_log,bins.(fn).centers_log] = ...
      SelectHistEdges(hist_params.nbins,bins.(fn).lims,'logabs');
  else
    [bins.(fn).edges_log,bins.(fn).centers_log] = ...
      SelectHistEdges(hist_params.nbins,bins.(fn).lims,'log');
  end
  
end

%% save results

savename = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframebinsfilestr);
save(savename,'bins');