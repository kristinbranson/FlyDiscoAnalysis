function FlyBowlComputePerFrameStats(expdir,varargin)

[analysis_protocol,settingsdir,datalocparamsfilestr,visible] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'visible','off'); %#ok<NASGU>

%% load this experiment
fprintf('Initializing trx...\n');

trx = Trx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir,...
  'datalocparamsfilestr',datalocparamsfilestr);

fprintf('Loading trajectories for %s...\n',expdir);

trx.AddExpDir(expdir);

%% read the stats params

statsperframefeaturesfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.statsperframefeaturesfilestr);
stats_perframefeatures = ReadStatsPerFrameFeatures(statsperframefeaturesfile);
statframeconditionsfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.statframeconditionfilestr);
frameconditiondict = ReadParams(statframeconditionsfile);
statflyconditionsfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.statflyconditionfilestr);
flyconditiondict = ReadParams(statflyconditionsfile);
statsparamsfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.statsparamsfilestr);
stats_params = ReadParams(statsparamsfile);

%% compute stats per fly

nflies = trx.nflies;
nprctiles = numel(stats_params.prctiles_compute);

statstxtsavename = fullfile(expdir,trx.dataloc_params.statsperframetxtfilestr);
statsmatsavename = fullfile(expdir,trx.dataloc_params.statsperframematfilestr);

% open the text file for writing
% statstxtfid = fopen(statstxtsavename,'w');
% fclose(statstxtfid);

statsperfly = struct;
statsperexp = struct;

for i = 1:numel(stats_perframefeatures),
  
  % which per-frame feature
  field = stats_perframefeatures(i).field;

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
  fn = sprintf('%s_fly%s_frame%s',field,flyconditionname,frameconditionname);
  
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
    data = trx(fly).(field);
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
      'prctiles_compute',stats_params.prctiles_compute);
    
  end

  statsperfly.(fn) = statsperflycurr;
  
  % compute per-experiment stats
  statsperexp.(fn) = CombinePerFrameStats(statsperflycurr);
  
  nflies_analyzed = nnz(~isnan(statsperfly.(fn).Z));
  nfliestotal = numel(statsperfly.(fn).Z);
  fprintf('Analyzed %d / %d flies for condition %s, %s\n',nflies_analyzed,nfliestotal,frameconditionname,flyconditionname);
  
  % save to text file
  %SavePerFrameStatsTxtFile(statstxtsavename,fn,statsperflycurr,statsperexp.(fn));
  
end

SaveAllPerFrameStatsTxtFile(statstxtsavename,statsperfly,statsperexp);

% save to mat file
save(statsmatsavename,'statsperfly','statsperexp','frameconditiondict',...
  'flyconditiondict','stats_params');

%% read the histogram params

histperframefeaturesfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.histperframefeaturesfilestr);
hist_perframefeatures = ReadHistPerFrameFeatures(histperframefeaturesfile);
% histparamsfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.histparamsfilestr);
% hist_params = ReadParams(histparamsfile);
histperframebinsfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.histperframebinsfilestr);
load(histperframebinsfile,'bins');

%% compute per fly histograms


histtxtsavename = fullfile(expdir,trx.dataloc_params.histperframetxtfilestr);
histmatsavename = fullfile(expdir,trx.dataloc_params.histperframematfilestr);

% open the text file for writing
% histtxtfid = fopen(histtxtsavename,'w');
% fclose(histtxtfid);

histperfly = struct;
histperexp = struct;

for i = 1:numel(hist_perframefeatures),
  
  % parameters for this histogram
  field = hist_perframefeatures(i).field;
  frameconditionname = hist_perframefeatures(i).framecondition;
  if ~strcmpi(frameconditionname,'any'),
    frameconditionparams = frameconditiondict.(frameconditionname);
  else
    frameconditionparams = {};
  end
  flyconditionname = hist_perframefeatures(i).flycondition;
  if ~strcmpi(flyconditionname,'any'),
    flyconditionparams = flyconditiondict.(flyconditionname);
  else
    flyconditionparams = {};
  end
  minZ = hist_perframefeatures(i).minZ;
  
  fn = sprintf('%s_fly%s_frame%s',field,flyconditionname,frameconditionname);

  fprintf('Computing histograms for %s...\n',fn);
  
  edges_linear = [-inf,bins.(field).edges_linear,inf];
  edges_log = [-inf,bins.(field).edges_log,inf];
  edges_linear(end-1) = edges_linear(end-1)*(1+1e-3);
  edges_log(end-1) = edges_log(end-1)*(1+1e-3);
  nbins_linear = numel(bins.(field).edges_linear)-1;
  nbins_log = numel(bins.(field).edges_log)-1;
  
  histperflycurr = struct('Z',nan(1,nflies),...
    'frac_linear',nan(nbins_linear,nflies),...
    'frac_log',nan(nbins_log,nflies),...
    'fracless_linear',nan(1,nflies),...
    'fracmore_linear',nan(1,nflies),...
    'fracless_log',nan(1,nflies),...
    'fracmore_log',nan(1,nflies),...
    'fracframesanalyzed',zeros(1,nflies));

  for fly = 1:nflies,
    
    % check that the fly matches the conditions in flyconditionparams
    if ~strcmpi(flyconditionname,'any') && ...
        ~FlyConditionCheck(trx,fly,flyconditionparams),
      %fprintf('Skipping fly %d for condition %s\n',fly,flyconditionname);
      continue;
    end
    
    data = trx(fly).(field);
    n = numel(data);
    
    % choose frames that match the conditions in frameconditionparams
    if strcmpi(frameconditionname,'any'),
      doanalyze = true(1,n);
    else
      doanalyze = FrameConditionCheck(trx,fly,n,frameconditionparams);
    end

    % make sure there are enough frames
    Z = nnz(doanalyze);
    if Z < minZ,
      %fprintf('Skipping fly %d for condition %s; not enough frames of data\n',fly,frameconditionname);
      continue;
    end
    
    % histogram
    [histperflycurr.frac_linear(:,fly),...
      histperflycurr.fracless_linear(fly),...
      histperflycurr.fracmore_linear(fly),...
      histperflycurr.Z(fly),...
      histperflycurr.fracframesanalyzed(fly)] = ...
      HistPerFrameData(data,doanalyze,edges_linear);

    [histperflycurr.frac_log(:,fly),...
      histperflycurr.fracless_log(fly),...
      histperflycurr.fracmore_log(fly)] = ...
      HistPerFrameData(data,doanalyze,edges_log);  
    
  end

  histperfly.(fn) = histperflycurr;
  % compute per-experiment hist
  histperexp.(fn) = CombinePerFrameHists(histperflycurr);
    
  nflies_analyzed = nnz(~isnan(histperfly.(fn).Z));
  nfliestotal = numel(histperfly.(fn).Z);
  fprintf('Analyzed %d / %d flies for condition %s, %s\n',nflies_analyzed,nfliestotal,frameconditionname,flyconditionname);

  
  % save to text file
  %SavePerFrameHistTxtFile(histtxtsavename,fn,histperflycurr,histperexp.(fn));

end

save(histmatsavename,'histperfly','histperexp','bins','frameconditiondict',...
  'flyconditiondict');
SaveAllPerFrameHistTxtFile(histtxtsavename,histperfly,histperexp);

% %% create the plot directory if it does not exist
% figdir = fullfile(expdir,trx.dataloc_params.figdir);
% if ~exist(figdir,'file'),
%   [status,msg,~] = mkdir(figdir);
%   if ~status,
%     error('Could not create the figure directory %s:\n%s',figdir,msg);
%   end
% end
% 
% %% read plotting parameters
% 
% histplotparamsfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.histplotparamsfilestr);
% hist_plot_params = ReadParams(histplotparamsfile);
% [~,expname] = fileparts(expdir);
% 
% %% plot stuff
% % TODO: plot stats of controls for the previous X weeks behind these
% 
% hist_fields = unique({hist_perframefeatures.field});
% for i = 1:numel(hist_fields),
%   
%   field = hist_fields{i};
% 
%   hfig = PlotPerFrameHists(field,hist_perframefeatures,...
%     histperexp,histperfly,...
%     bins.(field),hist_plot_params,expname,...
%     'visible',visible);
%   drawnow;
%   savename = sprintf('hist_%s.png',hist_fields{i});
%   if exist(savename,'file'),
%     delete(savename);
%   end
%   save2png(fullfile(figdir,savename),hfig);
%   
% end
% 
% close all;

