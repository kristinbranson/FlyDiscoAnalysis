function [perflystats,expstats] = ComputeHackHitCategoryStatistics(expdir,varargin)

[analysis_protocol,settingsdir,datalocparamsfilestr] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt');

%% read in the data locations
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

%% read in parameters
stats_params = ReadParams(fullfile(settingsdir,analysis_protocol,dataloc_params.hackhitcategoryparamsfilestr));

%% create trx variable

trx = FBATrx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir,...
  'datalocparamsfilestr',datalocparamsfilestr);

fprintf('Loading trajectories for %s...\n',expdir);

trx.AddExpDir(expdir);

%% initialize structures

statsperfly = struct;
statsperexp = struct;

%% compute fraction of time "jumping"

statsperfly.isjumping = struct;
statsperflycurr = struct;
for i = 1:trx.nflies,
  if trx(i).nframes < stats_params.min_nframes_jump,
    continue;
  end
  isjumping = double(trx(i).velmag >= stats_params.min_velmag_jump);
  doanalyze = true(size(isjumping));
  [statsperflycurr.Z(i),statsperflycurr.mean(i),...
    statsperflycurr.std(i),statsperflycurr.prctiles(:,i)] = ...
    ComputePerFrameStats(isjumping,doanalyze,...
    'prctiles_compute',[]);
  statsperflycurr.fracframesanalyzed(i) = nnz(doanalyze) / numel(doanalyze);
end
statsperexp.isjumping = CombinePerFrameStats(statsperflycurr);
statsperfly.isjumping = statsperflycurr;

%% compute mean speed while moving

statsperfly.meanspeedmoving = struct;
statsperflycurr = struct;
for i = 1:trx.nflies,
  if trx(i).nframes < stats_params.min_nframes_meanspeedmoving,
    continue;
  end
  velmag = trx(i).velmag;
  doanalyze = velmag >= stats_params.min_velmag_moving & ...
    velmag  <= stats_params.max_velmag_moving;
  [statsperflycurr.Z(i),statsperflycurr.mean(i),...
    statsperflycurr.std(i),statsperflycurr.prctiles(:,i)] = ...
    ComputePerFrameStats(velmag,doanalyze);
  statsperflycurr.fracframesanalyzed(i) = 1;
end
statsperexp.meanspeedmoving = CombinePerFrameStats(statsperflycurr);
statsperfly.meanspeedmoving = statsperflycurr;

%% save results to file

% SaveAllPerFrameStatsTxtFile(statstxtsavename,statsperfly,statsperexp);

% save to mat file
statsmatsavename = fullfile(expdir,dataloc_params.hackhitstatsmatfilestr);
save(statsmatsavename,'statsperfly','statsperexp','stats_params');