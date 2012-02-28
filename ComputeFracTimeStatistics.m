function [statsperfly,statsperexp] = ComputeFracTimeStatistics(expdir,varargin)

[analysis_protocol,settingsdir,datalocparamsfilestr] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt');

%% read parameters

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

stats_perframefeatures = ReadStatsPerFrameFeatures(fullfile(settingsdir,analysis_protocol,dataloc_params.behaviordetectionparamsfilestr));
statframeconditionsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statframeconditionfilestr);
frameconditiondict = ReadParams(statframeconditionsfile);
statflyconditionsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statflyconditionfilestr);
flyconditiondict = ReadParams(statflyconditionsfile);

%% create trx variable

trx = Trx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir,...
    'datalocparamsfilestr',datalocparamsfilestr);

fprintf('Loading trajectories for %s...\n',expdir);

trx.AddExpDir(expdir,'dooverwrite',false);

%% compute stats per fly

nflies = trx.nflies;

statstxtsavename = fullfile(expdir,trx.dataloc_params.fractimestatstxtfilestr);
statsmatsavename = fullfile(expdir,trx.dataloc_params.fractimestatsmatfilestr);

statsperfly = struct;
statsperexp = struct;
nprctiles = 0;

for i = 1:numel(stats_perframefeatures),
  
  % which behavior
  field = stats_perframefeatures(i).field;
  
  % check for mat file
  matfilestr = [field,'_labels.mat'];
  matfile = fullfile(expdir,matfilestr);
  if ~exist(matfile,'file'),
    warning('File %s does not exist',matfile);
    % TODO: add dummy values
    continue;
  end
  
  % load in labels
  labelidx = LoadLabelsFromFile(matfile,trx.nframes,trx.firstframes);
  
  % which frames to analyze
  frameconditionname = stats_perframefeatures(i).framecondition;
  if ~strcmpi(frameconditionname,'any'),
    frameconditionparams = frameconditiondict.(frameconditionname);
  else
    frameconditionparams = {};
  end

  % which flies to analyze
  flyconditionname = stats_perframefeatures(i).flycondition;
  if ~strcmpi(flyconditionname,'any'),
    flyconditionparams = flyconditiondict.(flyconditionname);
  else
    flyconditionparams = {};
  end
  
  % minimum number of frames to use this trajectory as a fly
  minZ = stats_perframefeatures(i).minZ;

  % unique identifier for this computation, composed of field + conditions
  fn = sprintf('fractime_%s_fly%s_frame%s',field,flyconditionname,frameconditionname);
  
  fprintf('Computing stats for %s...\n',fn);
  
  % initialize per-fly stats
  statsperflycurr = struct('Z',nan(1,nflies),...
    'mean',nan(1,nflies),...
    'std',nan(1,nflies),...
    'prctiles',nan(nprctiles,nflies),...
    'fracframesanalyzed',zeros(1,nflies));

  for fly = 1:nflies,
  
    % check that the fly matches the conditions in flyconditionparams
    if ~strcmpi(flyconditionname,'any') && ...
        ~FlyConditionCheck(trx,fly,flyconditionparams),
      %fprintf('Skipping fly %d for condition %s\n',fly,flyconditionname);
      continue;
    end

    % current field data
    data = double(labelidx{fly});
    n = numel(data);
    
    % choose frames that match the conditions in frameconditionparams
    if strcmpi(frameconditionname,'any'),
      doanalyze = true(1,n);
    else
      doanalyze = FrameConditionCheck(trx,fly,n,frameconditionparams);
    end
 
    % skip this trajectory if there aren't enough frames of data
    if nnz(doanalyze) < minZ,
      %fprintf('Skipping fly %d for condition %s; not enough frames of data\n',fly,frameconditionname);
      continue;
    end

    % we will weight this trajectory by the fraction of frames analyzed
    statsperflycurr.fracframesanalyzed(fly) = nnz(doanalyze) / n;
    
    % compute stats
    [statsperflycurr.Z(fly),statsperflycurr.mean(fly),...
      statsperflycurr.std(fly),statsperflycurr.prctiles(:,fly)] = ...
      ComputePerFrameStats(data,doanalyze,...
      'prctiles_compute',[]);

  end
  
  statsperfly.(fn) = statsperflycurr;
  
  % compute per-experiment stats
  statsperexp.(fn) = CombinePerFrameStats(statsperflycurr);
  
  nflies_analyzed = nnz(~isnan(statsperfly.(fn).Z));
  nfliestotal = numel(statsperfly.(fn).Z);
  fprintf('Analyzed %d / %d flies for condition %s, %s\n',nflies_analyzed,nfliestotal,frameconditionname,flyconditionname);
    
end

SaveAllPerFrameStatsTxtFile(statstxtsavename,statsperfly,statsperexp);

% save to mat file
save(statsmatsavename,'statsperfly','statsperexp','frameconditiondict',...
  'flyconditiondict');
