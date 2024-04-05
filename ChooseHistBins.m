%% set up path
modpath
%%
settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
realrootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
%rootdir = '/groups/branson/bransonlab/projects/olympiad/HackHitData';
%rootdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlData';
rootdir = realrootdir;
rootdir_fixed = '';

analysis_protocol = 'current';
datalocparamsfilestr = 'dataloc_params.txt';
params = {'analysis_protocol',analysis_protocol,...
  'datalocparamsfilestr',datalocparamsfilestr,...
  'settingsdir',settingsdir};

expdirfile = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/JAABA_guff/perframe/expdirs_jaabadetect20130606.txt';
%expdirfile = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/experiments/datamanagement/expdirs_primary20130306.txt';

behaviornames = {
  'attemptedcopulation'
  'backup'
  'bodyturn'
  'chase'
  'copulation'
  'crabwalkall'
  'crabwalkextreme'
  'jump'
  'pivotcenter'
  'pivottail'
  'righting'
  'stop'
  'touch'
  'wingextension'
  'winggrooming'
  'wingflick'
  'walk'
  };


% behaviornames = {
%   'walking'
%   'stopped'
%   'jumping'
%   'righting'
%   'chasing'
%   'backup'
%   'crabwalk'
%   'pivottail'
%   'pivothead'
%   'touch'
%   'winggrooming'
%   'copulation'
%   'attemptedcopulation'
%   };
%behaviornames = {'chase','walk','stop','jump'};

fps = 30.3445;

%% parameters

datalocparamsfile = fullfile(settingsdir,analysis_protocol,'dataloc_params.txt');
dataloc_params = ReadParams(datalocparamsfile);

%% choose experiments to histogram

moviefilestr = dataloc_params.moviefilestr;
trxfilestr = dataloc_params.ctraxfilestr;
annfilestr = dataloc_params.annfilestr;
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

% experiment_params = struct;
% % root data dir
% experiment_params.rootdir = rootdir;
% % what dates should we analyze
% %experiment_params.daterange = {'20110201T000000','20111221T000000'};
% % remove failures
% experiment_params.checkflags = false;
% % remove missing data
% experiment_params.removemissingdata = false;
% % azanchir experiments
% experiment_params.screen_type = 'non_olympiad_azanchir*';
% 
% % for specifying galits data
% experiment_params.line_name = 'EXT_CantonS_*';
% experiment_params.effector = 'NoEffector_0_9999';
% experiment_params.handling_protocol = ...
%   {'HP_flybowl_v007p1.xls','HP_flybowl_v007p2.xls','HP_flybowl_v007p3.xls',...
%   'HP_flybowl_v007p6.xls','HP_flybowl_v007p7.xls'};
% experiment_params.protocol = 'EP_flybowl_v011p3.xls';
% 
% % how to weight various things
weight_order = {'genotype','date','rig','bowl'};
% % required files
% %subreadfiles = {'chase_labels.mat','stop_labels.mat','walk_labels.mat','jump_labels.mat'};
% subreadfiles = {'perframe/dnose2ell_angle_min20to20.mat'};
% % what lines
% %experiment_params.linename = '';
% % 
ngal4 = 100;
ncontrols = 100;
% tmpparams = struct2paramscell(experiment_params);
% 
% %[~,~,~,experiments_all] = getExperimentDirs(tmpparams{:});
% experiments_all = SAGEListBowlExperiments(tmpparams{:});

expdirs_all = importdata(expdirfile);
clear experiments_all;
for i = 1:numel(expdirs_all),
  metadatacurr = parseExpDir(expdirs_all{i},true);
  metadatacurr.line = metadatacurr.line_name;
  experiments_all(i) = metadatacurr;
end


% baddata = false(1,numel(experiments_all));
% for i = 1:numel(experiments_all),
%   if ~exist(experiments_all(i).file_system_path,'dir'),
%     baddata(i) = true;
%     continue;
%   end
%   [~,experiment_name] = myfileparts(experiments_all(i).file_system_path);
% %   annfile = fullfile(rootdir,experiment_name,annfilestr);
% %   if ~exist(annfile,'file'),
% %     baddata(i) = true;
% %   end
%   for j = 1:numel(subreadfiles),
%     if ~exist(fullfile(rootdir,experiment_name,subreadfiles{j}),'file'),
%       baddata(i) = true;
%     end
%   end
% end
% experiments_all(baddata) = [];

% 
% baddata = false(1,numel(experiments_all));
% mindatenum = datenum('20120712','yyyymmdd');
% for i = 1:numel(experiments_all),
%   fprintf('%d: %s\n',i,experiments_all(i).experiment_name);
%   tmp = dir(fullfile(experiments_all(i).file_system_path,'registered_trx.mat'));
%   if tmp.datenum < mindatenum,
%     baddata(i) = true;
%     fprintf('dt not fixed for %s\n',experiments_all(i).experiment_name);
%   end
%   
% %   tmp = who('-file',fullfile(experiments_all(i).file_system_path,'registered_trx.mat'));
% %   if ~ismember('dtfixed',tmp),
% %     fprintf('dt not fixed for %s\n',experiments_all(i).experiment_name);
% %     baddata(i) = true;
% %   end
% end
%     

%experiments_all = rmfield(experiments_all,'date');

[experiments,ngal4_chosen,ncontrols_chosen] = choose_expdirs(experiments_all,ngal4,ncontrols,'weight_order',weight_order);

expdirs = {experiments.file_system_path};
nexpdirs = numel(expdirs);



%% load data

% obj = FBATrx('settingsdir',settingsdir,'analysis_protocol',analysis_protocol);
% for i = 1:numel(expdirs),
%   obj.AddExpDir(expdirs{i});
% end

%% read parameters

% perframefnsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.perframefnsfilestr);
% perframefns = importdata(perframefnsfile);
prctile_store = 10;

histparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histparamsfilestr);

histperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframefeaturesfilestr);
hist_perframefeatures = ReadHistPerFrameFeatures2(histperframefeaturesfile);
hist_params = ReadParams(histparamsfile);

statframeconditionsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statframeconditionfilestr);
frameconditiondict = ReadParams(statframeconditionsfile);
statflyconditionsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statflyconditionfilestr);
flyconditiondict = ReadParams(statflyconditionsfile);



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
      m = m{1};
      s = regexp(frameconditiondict.(m){1},'^[^_]+_(.+)$','tokens','once');
      if isempty(s),
        error('Could not parse %s',m);
      end
      filename = fullfile(expdir,[s{1},'.mat']);
      if exist(filename,'file'),
        labeldata = load(filename);
        datacurr = (cell2mat(labeldata.t1s)-cell2mat(labeldata.t0s))/fps;
      else
        s{1} = regexprep(s{1},'label','score');
        filename = fullfile(expdir,[s{1},'.mat']);
        if ~exist(filename,'file'),
          error('Could not find score or label file for %s',m);
        end
        scoredata = load(filename);
        if isfield(scoredata.allScores,'postprocessed'),
          datacurr = [];
          for j = 1:numel(scoredata.allScores.postprocessed),
            scoredata.allScores.postprocessed{j}(isnan(scoredata.allScores.postprocessed{j})) = 0;
            [t0s,t1s] = get_interval_ends(scoredata.allScores.postprocessed{j});
            datacurr = [datacurr,(t1s-t0s)/fps]; %#ok<AGROW>
          end
        else
          datacurr = (cell2mat(scoredata.allScores.t1s)-cell2mat(scoredata.allScores.t0s))/fps;
        end
      end
    else
      filename = fullfile(expdir,dataloc_params.perframedir,[fn,'.mat']);
      if ~exist(filename,'file'), 
        warning('Could not find file %s, skipping',filename);
        continue;
      end
      % get per-frame data for these flies
      tmp = load(filename,'data');
      if ~iscell(tmp.data),
        warning('Data in %s is not a cell, skipping.',filename);
        continue;
      end
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