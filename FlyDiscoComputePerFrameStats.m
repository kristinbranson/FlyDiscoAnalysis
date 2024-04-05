function FlyDiscoComputePerFrameStats(expdir,varargin)
% debugplot,n will create figures that display perframe data and the time period of the first n statistics
% defined in stats_perframefeatures.txt with the norm periods, if any, plotted in
% the subplot below 
QUICKDEBUG = false;

if QUICKDEBUG,
  global trx; %#ok<UNRCH,TLEV>
  warning('using global variable trx');
end

version = '0.1';

special_cases = {'fractime','duration','boutfreq'};

[analysis_protocol,settingsdir,datalocparamsfilestr,...
  dorecompute,debug,verbose,docomputehists,debugplot] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'dorecompute',true,...
  'debug',false,...
  'verbose',true,...
  'docomputehists',true,...
  'debugplot',0); 

if ischar(docomputehists),
  docomputehists = str2double(docomputehists);
end
if ischar(verbose),
  verbose = str2double(verbose);
end
if ischar(debug),
  debug = str2double(debug);
end
if ischar(dorecompute),
  dorecompute = str2double(dorecompute);
end
if ischar(debugplot),
  debugplot = str2double(debugplot);
end

if verbose,
  fprintf('Analysis protocol: %s\n',analysis_protocol);
end

%% load this experiment
if verbose,
  fprintf('Initializing trx...\n');
end

if QUICKDEBUG && ~isempty(trx), %#ok<NODEF>
  warning('using global variable trx');
else
trx = FBATrx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir,...
  'datalocparamsfilestr',datalocparamsfilestr,...
  'maxdatacached',2^30);

if verbose,
  fprintf('Loading trajectories for %s...\n',expdir);
end

trx.AddExpDir(expdir,'dooverwrite',false,'openmovie',false);
end

% % load labels
% labelfiles = dir(fullfile(trx.expdirs{1},'*_labels.mat'));
% labelfiles = {labelfiles.name};
% requiredfns = {'flies','names','off','t0s','t1s'};
% for i = 1:numel(labelfiles),
%   tmp = who('-file',fullfile(trx.expdirs{1},labelfiles{i}));
%   if ~isempty(setdiff(requiredfns,tmp)),
%     continue;
%   end
%   trx.LoadLabelsFromFile(labelfiles{i});
% end

%% log file

if isfield(trx.dataloc_params,'perframestats_logfilestr') && ~debug,
  logfile = fullfile(expdir,trx.dataloc_params.perframestats_logfilestr);
  logfid = fopen(logfile,'a');
  if logfid < 1,
    warning('Could not open log file %s\n',logfile);
    logfid = 1;
  end
else
  logfid = 1;
end

timestamp = datestr(now,'yyyymmddTHHMMSS');
real_analysis_protocol = GetRealAnalysisProtocol(analysis_protocol,settingsdir);

fprintf(logfid,'\n\n***\nRunning FlyDiscoComputePerFrameStats version %s analysis_protocol %s (linked to %s) at %s\n',version,analysis_protocol,real_analysis_protocol,timestamp);


%% read the stats params

statsperframefeaturesfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.statsperframefeaturesfilestr);
stats_perframefeatures = ReadStatsPerFrameFeatures2(statsperframefeaturesfile);
statframeconditionsfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.statframeconditionfilestr);
frameconditiondict = ReadParams(statframeconditionsfile);
statflyconditionsfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.statflyconditionfilestr);
flyconditiondict = ReadParams(statflyconditionsfile);
statsparamsfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.statsparamsfilestr);
stats_params = ReadParams(statsparamsfile);

%% compute stats per fly

nflies = trx.nflies;

statstxtsavename = fullfile(expdir,trx.dataloc_params.statsperframetxtfilestr);
statsmatsavename = fullfile(expdir,trx.dataloc_params.statsperframematfilestr);

if ~dorecompute && exist(statstxtsavename,'file') && exist(statsmatsavename,'file'),
  load(statsmatsavename); %#ok<LOAD>
else

% open the text file for writing
% statstxtfid = fopen(statstxtsavename,'w');
% fclose(statstxtfid);


% delete data we should recompute
if exist(statstxtsavename,'file'),
  try
    fprintf('Deleting file %s...\n',statstxtsavename);
    delete(statstxtsavename);
  catch
    warning('Could not delete file %s',statstxtsavename);
  end
end
if exist(statsmatsavename,'file'),
  try
    fprintf('Deleting file %s...\n',statsmatsavename);
    delete(statsmatsavename);
  catch
    warning('Could not delete file %s',statsmatsavename);
  end
end
perframefiles = dir(fullfile(expdir,trx.dataloc_params.perframefnsfilestr,'stim*.mat'));
perframefiles = {perframefiles.name};
perframefiles(cellfun(@isempty,regexp(perframefiles,'^stim((on)|(off))(_\d+)+\.mat','once'))) = [];
for i = 1:numel(perframefiles),
  try
    pff = fullfile(expdir,trx.dataloc_params.perframefnsfilestr,perframefiles{i});
    fprintf('Deleting file %s...\n',pff);
    delete(pff);
  catch
    warning('Could not delete file %s',pff);
  end
end

statsperfly = struct;
statsperexp = struct;

allfns = cell(1,numel(stats_perframefeatures));

if debugplot > 0,
  hfig = figure;
  hax = createsubplots(2,debugplot,.05,hfig);
  hax = reshape(hax,[2,debugplot]);
  for i = 1:numel(hax),
    hold(hax(i),'on');
  end
end  

for i = 1:numel(stats_perframefeatures),
  
  [fn,fullfn,field,flyconditionname,frameconditionname,frameconditionparams,...
    flyconditionparams,minZboth,minZfly,...
    isnorm,norm_flyconditionname,norm_frameconditionname,...
    norm_frameconditionparams,norm_flyconditionparams,norm_minZboth,norm_minZfly] = ...
    ParseStatParams(trx,stats_perframefeatures(i),allfns(1:i-1),special_cases,stats_params,...
    frameconditiondict,flyconditiondict);
  allfns{i} = fn;
  
  if verbose,
    fprintf('Computing stats for %s...\n',fn);
  end
  
  if i <= debugplot,
    haxcurr = hax(1,i);
    haxcurrnorm = hax(2,i);
  else
    haxcurr = [];
    haxcurrnorm = [];
  end
  
  [statsperflycurr,alldata,alldatanorm] = ComputeStat(trx,field,flyconditionname,frameconditionname,flyconditionparams,frameconditionparams,stats_params,special_cases,haxcurr);
  if ~isempty(haxcurr),
    title(haxcurr,sprintf('%s fly=%s frame=%s',field,flyconditionname,frameconditionname),'interpreter','none');
  end

  if isnorm,
    [norm_statsperflycurr,norm_alldata,norm_alldatanorm] = ComputeStat(trx,field,norm_flyconditionname,norm_frameconditionname,norm_flyconditionparams,norm_frameconditionparams,stats_params,special_cases,haxcurrnorm);
    statsperflycurr = ComputeNormStat(statsperflycurr,norm_statsperflycurr);
    statsperexpcurr = CombinePerFrameStats2(statsperflycurr,minZboth,minZfly,field,alldata,alldatanorm,norm_minZboth,norm_minZfly,norm_alldata,norm_alldatanorm);
    if ~isempty(haxcurrnorm),
      title(haxcurrnorm,sprintf('%s fly=%s frame=%s',field,norm_flyconditionname,norm_frameconditionname),'interpreter','none');
    end
  else
    statsperexpcurr = CombinePerFrameStats2(statsperflycurr,minZboth,minZfly,field,alldata,alldatanorm);
  end
  statsperflycurr.name = fullfn;
  statsperexpcurr.name = fullfn;
    
  statsperfly.(fn) = statsperflycurr;
  statsperexp.(fn) = statsperexpcurr;
  
  %nflies_analyzed = nnz(~isnan(statsperfly.(fn).Z));
  % fracframesanalyzed corresponds to doanalyze_fly
  if verbose,
    
    if isnorm,
      
      idx = statsperfly.(fn).raw_Z > 0;
      raw_nflies_fly = sum(statsperfly.(fn).raw_fracframesanalyzed(idx));
      raw_nframes_per_fly = sum(statsperfly.(fn).raw_Z(idx)) / sum(raw_nflies_fly);
      raw_fracframes_per_fly = sum(statsperfly.(fn).raw_Z(idx)) / sum(statsperfly.(fn).raw_Zfly(idx));
      idx = statsperfly.(fn).norm_Z > 0;
      norm_nflies_fly = sum(statsperfly.(fn).norm_fracframesanalyzed(idx));
      norm_nframes_per_fly = sum(statsperfly.(fn).norm_Z(idx)) / sum(norm_nflies_fly);
      norm_fracframes_per_fly = sum(statsperfly.(fn).norm_Z(idx)) / sum(statsperfly.(fn).norm_Zfly(idx));
      if strcmp(field,'fractime'),
        fprintf('%.1f flies satisfy conditions %s, %.1f flies satisfy conditions %s, on average %.1f normed by %.1f frames per fly\n',raw_nflies_fly,flyconditionname,norm_nflies_fly,norm_flyconditionname,raw_nframes_per_fly,norm_nframes_per_fly);
      else
        fprintf('%.1f flies satisfy raw  conditions %s, on average %.1f frames (%.1f%%) analyzed per fly also satisfy %s\n',raw_nflies_fly,flyconditionname,raw_nframes_per_fly,raw_fracframes_per_fly*100,frameconditionname);
        fprintf('%.1f flies satisfy norm conditions %s, on average %.1f frames (%.1f%%) analyzed per fly also satisfy %s\n',norm_nflies_fly,norm_flyconditionname,norm_nframes_per_fly,norm_fracframes_per_fly*100,norm_frameconditionname);
      end
      
    else
      idx = statsperfly.(fn).Z > 0;
      nflies_fly = sum(statsperfly.(fn).fracframesanalyzed(idx));
      nframes_per_fly = sum(statsperfly.(fn).Z(idx)) / sum(nflies_fly);
      fracframes_per_fly = sum(statsperfly.(fn).Z(idx)) / sum(statsperfly.(fn).Zfly(idx));
      if strcmp(field,'fractime'),
        fprintf('%.1f flies satisfy conditions %s, on average %.1f frames analyzed per fly\n',nflies_fly,flyconditionname,nframes_per_fly);
      else
        fprintf('%.1f flies satisfy conditions %s, on average %.1f frames (%.1f%%) analyzed per fly also satisfy %s\n',nflies_fly,flyconditionname,nframes_per_fly,fracframes_per_fly*100,frameconditionname);
      end
    end
  end
  % save to text file
  %SavePerFrameStatsTxtFile(statstxtsavename,fn,statsperflycurr,statsperexp.(fn));
  
end

if ~debug,
  SaveAllPerFrameStatsTxtFile(statstxtsavename,statsperfly,statsperexp);
  
  % save to mat file
  if exist(statsmatsavename,'file'),
    try
      delete(statsmatsavename);
    catch
      warning('Could not remove file %s before saving to it.',statsmatsavename);
    end
  end
  save(statsmatsavename,'statsperfly','statsperexp','frameconditiondict',...
    'flyconditiondict','stats_params',...
    'version','timestamp','analysis_protocol','real_analysis_protocol');
end
end

%% hists

if docomputehists,


%% read the histogram params


histperframefeaturesfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.histperframefeaturesfilestr);
hist_perframefeatures = ReadStatsPerFrameFeatures2(histperframefeaturesfile);
% histparamsfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.histparamsfilestr);
% hist_params = ReadParams(histparamsfile);
histperframebinsfile = fullfile(trx.settingsdir,trx.analysis_protocol,trx.dataloc_params.histperframebinsfilestr);
load(histperframebinsfile,'bins');

%% compute per fly histograms


histtxtsavename = fullfile(expdir,trx.dataloc_params.histperframetxtfilestr);
histmatsavename = fullfile(expdir,trx.dataloc_params.histperframematfilestr);

if ~dorecompute && exist(histtxtsavename,'file') && exist(histmatsavename,'file'),
  load(histmatsavename); %#ok<LOAD>
else

% open the text file for writing
% histtxtfid = fopen(histtxtsavename,'w');
% fclose(histtxtfid);

histperfly = struct;
histperexp = struct;

nbins_linear_default = [];
nbins_log_default = [];

allhistfns = cell(1,numel(hist_perframefeatures));

for i = 1:numel(hist_perframefeatures),
  
  [fn,fullfn,field,flyconditionname,frameconditionname,frameconditionparams,...
    flyconditionparams,minZboth,minZfly,isnorm] = ...
    ParseStatParams(trx,hist_perframefeatures(i),allhistfns(1:i-1),special_cases,stats_params,...
    frameconditiondict,flyconditiondict);
  assert(~isnorm,'Normalization not implemented for histogram statistics');
  allhistfns{i} = fullfn;
  

  if verbose,
    fprintf('Computing histograms for %s...\n',fn);
  end
  
  if ismember(field,special_cases),
    binfn = [field,'_',frameconditionname];
%     tmp = frameconditionparams(1:2:end);
%     m = regexp(tmp,'^[^_]+_(.+)_labels$','once','tokens');
%     tmp1 = cellfun(@isempty,m);
%     m(tmp1) = regexp(tmp(tmp1),'^[^_]+_labels_(.+)$','once','tokens');
%     tmp = find(~cellfun(@isempty,m),1);
%     binfn = [field,'_',m{tmp}{1}];
  else
    binfn = field;
  end
  
  if ~isfield(bins,binfn),
    warning('Field %s missing from histogram bins',binfn);
    if isempty(nbins_linear_default),
      tmpfns = fieldnames(bins);
      nbins_linear_all = nan(1,numel(tmpfns));
      nbins_log_all = nan(1,numel(tmpfns));
      for tmpi = 1:numel(tmpfns),
        tmpfn = tmpfns{tmpi};
        nbins_linear_all(tmpi) = numel(bins.(tmpfn).centers_linear);
        nbins_log_all(tmpi) = numel(bins.(tmpfn).centers_log);
        nbins_linear_default = mode(nbins_linear_all);
        nbins_log_default = mode(nbins_log_all);
      end
    end
    nbins_linear = nbins_linear_default;
    nbins_log = nbins_log_default;
    
    alldata = [];
  
    for fly = 1:nflies,
      
      % current field data
      if ~ismember(field,special_cases),
        data = trx(fly).(field);
        n = numel(data);
      else
        n = trx(fly).nframes;
      end
      
      % check that the fly matches the conditions in flyconditionparams
      if strcmpi(flyconditionname,'any'),
        doanalyze_fly = true(1,n);
      else
        doanalyze_fly = FrameConditionCheck(trx,fly,n,flyconditionparams);
      end
      
      % choose frames that match the conditions in frameconditionparams
      if strcmpi(frameconditionname,'any'),
        doanalyze_frame = true(1,n);
      else
        doanalyze_frame = FrameConditionCheck(trx,fly,n,frameconditionparams);
      end
      
      doanalyze = doanalyze_fly & doanalyze_frame;
      
      % skip this trajectory if there aren't enough frames of data
      if nnz(doanalyze) < minZboth || nnz(doanalyze_fly) < minZfly,
        %fprintf('Skipping fly %d for condition %s; not enough frames of data\n',fly,frameconditionname);
        continue;
      end
      
      if strcmp(field,'duration'),
        [bout_starts,bout_ends] = ComputeBouts(doanalyze,doanalyze_fly);
        data = trx(fly).timestamps(min(bout_ends+1,n))-trx(fly).timestamps(bout_starts);
        alldata = [alldata,data]; %#ok<AGROW>
      else
        alldata = [alldata,data(doanalyze)]; %#ok<AGROW>
      end
    end

    lim = [min(alldata),max(alldata)];
    edges_linear = SelectHistEdges(nbins_linear,lim,'linear');
    if lim(1) < 0,
      edges_log = SelectHistEdges(nbins_log,lim,'logabs');
    else
      edges_log = SelectHistEdges(nbins_log,lim,'log');
    end
    edges_linear = [-inf,edges_linear,inf]; %#ok<AGROW>
    edges_log = [-inf,edges_log,inf]; %#ok<AGROW>
    edges_linear(end-1) = edges_linear(end-1)*(1+1e-3);
    edges_log(end-1) = edges_log(end-1)*(1+1e-3);
    
  else
    edges_linear = [-inf,bins.(binfn).edges_linear,inf];
    edges_log = [-inf,bins.(binfn).edges_log,inf];
    edges_linear(end-1) = edges_linear(end-1)*(1+1e-3);
    edges_log(end-1) = edges_log(end-1)*(1+1e-3);
    nbins_linear = numel(bins.(binfn).edges_linear)-1;
    nbins_log = numel(bins.(binfn).edges_log)-1;
  end
  
  histperflycurr = struct('Z',zeros(1,nflies),...
    'frac_linear',nan(nbins_linear,nflies),...
    'frac_log',nan(nbins_log,nflies),...
    'fracless_linear',nan(1,nflies),...
    'fracmore_linear',nan(1,nflies),...
    'fracless_log',nan(1,nflies),...
    'fracmore_log',nan(1,nflies),...
    'fracframesanalyzed',zeros(1,nflies),...
    'Zfly',zeros(1,nflies));

  %alldata = [];
%   
%   if isnorm,
%     norm_statsperflycurr = ComputeStat(trx,field,norm_flyconditionname,norm_frameconditionname,norm_flyconditionparams,norm_frameconditionparams,stats_params,special_cases);
%   end
%   
  for fly = 1:nflies,
    
    [data,doanalyze,doanalyze_fly,n] = ...
      ComputeDataFly(trx,fly,field,flyconditionname,frameconditionname,flyconditionparams,frameconditionparams,special_cases);
    
    % skip this trajectory if there aren't enough frames of data
    if nnz(doanalyze) < minZboth || nnz(doanalyze_fly) < minZfly,
      %fprintf('Skipping fly %d for condition %s; not enough frames of data\n',fly,frameconditionname);
      continue;
    end

    if strcmp(field,'duration'),
      [bout_starts,bout_ends] = ComputeBouts(doanalyze,doanalyze_fly);
      data = trx(fly).timestamps(min(bout_ends+1,n))-trx(fly).timestamps(bout_starts);
      %alldata = [alldata,data]; %#ok<AGROW>
      doanalyze_data = true(1,numel(data));
    else
      doanalyze_data = doanalyze;
      %alldata = [alldata,data(doanalyze)]; %#ok<AGROW>
    end
    
%     if isnorm,
%       data = data - norm_statsperflycurr.mean(fly);
%     end
    
    % histogram
    [histperflycurr.frac_linear(:,fly),...
      histperflycurr.fracless_linear(fly),...
      histperflycurr.fracmore_linear(fly),...
      histperflycurr.Z(fly)] = ...
      HistPerFrameData(data,doanalyze_data,edges_linear);

    [histperflycurr.frac_log(:,fly),...
      histperflycurr.fracless_log(fly),...
      histperflycurr.fracmore_log(fly)] = ...
      HistPerFrameData(data,doanalyze_data,edges_log);  
    
    % we will weight this trajectory by the fraction of frames analyzed
    histperflycurr.Zfly(fly) = nnz(doanalyze_fly);
    histperflycurr.fracframesanalyzed(fly) = nnz(doanalyze_fly) / n;
    histperflycurr.name = fullfn;
    
  end

  histperfly.(fn) = histperflycurr;
  % compute per-experiment hist
  histperexp.(fn) = CombinePerFrameHists2(histperflycurr,hist_perframefeatures(i).minZboth,hist_perframefeatures(i).minZfly);
  histperexp.(fn).name = fullfn;
  idx = histperfly.(fn).Z > 0;
  nflies_fly = sum(histperfly.(fn).fracframesanalyzed(idx));
  nframes_per_fly = sum(histperfly.(fn).Z(idx)) / sum(nflies_fly);
  fracframes_per_fly = sum(histperfly.(fn).Z(idx)) / sum(histperfly.(fn).Zfly(idx));
  if verbose,
    fprintf('Histogramming: %.1f flies satisfy conditions %s, on average %.1f frames (%.1f%%) analyzed per fly also satisfy %s\n',nflies_fly,flyconditionname,nframes_per_fly,fracframes_per_fly*100,frameconditionname);
  end
  
  % save to text file
  %SavePerFrameHistTxtFile(histtxtsavename,fn,histperflycurr,histperexp.(fn));

end

if ~debug,
  if exist(histmatsavename,'file'),
    try %#ok<TRYNC>
      delete(histmatsavename);
    end
  end
  save(histmatsavename,'histperfly','histperexp','bins','frameconditiondict',...
    'flyconditiondict','histperframefeaturesfile','hist_perframefeatures',...
    'version','timestamp','analysis_protocol','real_analysis_protocol');
  SaveAllPerFrameHistTxtFile(histtxtsavename,histperfly,histperexp);
end

end

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

end

%% close log

fprintf(logfid,'Finished running FlyDiscoComputePerFrameStats at %s.\n',datestr(now,'yyyymmddTHHMMSS'));

if logfid > 1,
  fclose(logfid);
end

function [fn,fullfn,field,flyconditionname,frameconditionname,frameconditionparams,...
    flyconditionparams,minZboth,minZfly,...
    isnorm,norm_flyconditionname,norm_frameconditionname,...
    norm_frameconditionparams,norm_flyconditionparams,norm_minZboth,norm_minZfly] = ...
  ParseStatParams(trx,param,allfns,special_cases,stats_params,frameconditiondict,flyconditiondict)

% which per-frame feature
field = param.field;

% which frames to analyze
frameconditionname = param.framecondition;
frameconditionparams = DecodeConditions(frameconditionname,frameconditiondict);

% which flies to analyze
flyconditionname = param.flycondition;
flyconditionparams = DecodeConditions(flyconditionname,flyconditiondict);

isnorm = ~isempty(param.norm_framecondition);
norm_flyconditionname = '';
norm_frameconditionname = '';
norm_frameconditionparams = [];
norm_flyconditionparams = [];
if isnorm,
  norm_frameconditionname = param.norm_framecondition;
  norm_frameconditionparams = DecodeConditions(norm_frameconditionname,frameconditiondict);
  norm_flyconditionname = param.norm_flycondition;
  norm_flyconditionparams = DecodeConditions(norm_flyconditionname,flyconditiondict);
end


% is this a multi-frame feature?
% TODO: this only works on features that depend on at most pairs of frames
[~,fly] = max(trx.nframes);
if ismember(field,special_cases) || numel(trx(fly).(field)) < trx.nframes(fly),
  % add in bounds on time between frames
  if isfield(stats_params,'min_dt'),
    flyconditionparams(end+1:end+2) = {'min_dt',num2str(stats_params.min_dt,10)};
  end
  if isfield(stats_params,'max_dt'),
    flyconditionparams(end+1:end+2) = {'max_dt',num2str(stats_params.max_dt,10)};
  end
  if isnorm,
    if isfield(stats_params,'min_dt'),
      norm_flyconditionparams(end+1:end+2) = {'min_dt',num2str(stats_params.min_dt,10)};
    end
    if isfield(stats_params,'max_dt'),
      norm_flyconditionparams(end+1:end+2) = {'max_dt',num2str(stats_params.max_dt,10)};
    end
  end
end

% minimum number of frames to use this trajectory as a fly
minZboth = param.minZboth;
minZfly = param.minZfly;

norm_minZboth = [];
norm_minZfly = [];
if isnorm,
  % minimum number of frames to use this trajectory as a fly
  norm_minZboth = param.norm_minZboth;
  norm_minZfly = param.norm_minZfly;
end

% unique identifier for this computation, composed of field + conditions
fn = sprintf('%s_fly%s_frame%s',field,flyconditionname,frameconditionname);
if isnorm,
  fn = sprintf('%s_norm_fly%s_frame%s',fn,norm_flyconditionname,norm_frameconditionname);
end
fullfn = fn;
if numel(fn) > 63,
  fn = fn(1:61);
  nprev = nnz(strcmp(allfns,fn));
  fn = [fn,num2str(nprev+1)]; 
end


function [statsperflycurr,alldata,alldatanorm] = ComputeStat(trx,field,flyconditionname,frameconditionname,flyconditionparams,frameconditionparams,stats_params,special_cases,haxcurr)

if ~exist('haxcurr','var'),
  haxcurr = [];
end

nflies = trx.nflies;
nprctiles = numel(stats_params.prctiles_compute);

% initialize per-fly stats
statsperflycurr = struct('Z',zeros(1,nflies),...
  'mean',nan(1,nflies),...
  'std',nan(1,nflies),...
  'prctiles',nan(nprctiles,nflies),...
  'fracframesanalyzed',zeros(1,nflies),...
  'Zfly',zeros(1,nflies));

if ismember(field,{'fractime','boutfreq'}),
  alldata = 0;
  alldatanorm = 0;
else
  alldata = [];
  alldatanorm = [];
end

hcurr = [];
[~,longestfly] = max(trx.nframes);
for fly = 1:nflies,
  
  if isempty(haxcurr) || fly~=longestfly,
    haxcurr1 = [];
  else
    haxcurr1 = haxcurr;
  end
  [data,doanalyze,doanalyze_fly,n,hcurr1] = ...
    ComputeDataFly(trx,fly,field,flyconditionname,frameconditionname,...
    flyconditionparams,frameconditionparams,special_cases,haxcurr1);
  
  if ~isempty(haxcurr1),
    hcurr = [hcurr,hcurr1]; %#ok<AGROW>
    plot(haxcurr1,trx(fly).firstframe:trx(fly).firstframe+numel(data)-1,data,'k.');
    plot(haxcurr1,trx(fly).firstframe - 1 + find(doanalyze&doanalyze_fly),data(doanalyze&doanalyze_fly),'.','Color',[0.4940    0.1840    0.5560]);
    ylabel(haxcurr1,field,'interpreter','none');
  end
  
  if strcmp(field,'fractime'),
    statsperflycurr.mean(fly) = nnz(doanalyze) / nnz(doanalyze_fly);
    % note that this is different!
    statsperflycurr.Z(fly) = nnz(doanalyze_fly);
    statsperflycurr.std(fly) = nan;
    statsperflycurr.prctiles(:,fly) = nan;
    alldata = alldata + nnz(doanalyze);
    alldatanorm = alldatanorm + nnz(doanalyze_fly);
    
  elseif strcmp(field,'boutfreq'),
    [bout_starts] = ComputeBouts(doanalyze,doanalyze_fly);
    nbouts = numel(bout_starts);
    statsperflycurr.mean(fly) = nbouts / nnz(doanalyze_fly);
    statsperflycurr.Z(fly) = nnz(doanalyze_fly);
    statsperflycurr.std(fly) = nan;
    statsperflycurr.prctiles(:,fly) = nan;
    alldata = alldata + nbouts;
    alldatanorm = alldatanorm + nnz(doanalyze_fly);
    
  elseif strcmp(field,'duration'),
    [bout_starts,bout_ends] = ComputeBouts(doanalyze,doanalyze_fly);
    boutdurs = trx(fly).timestamps(min(bout_ends+1,n))-trx(fly).timestamps(bout_starts);
    statsperflycurr.mean(fly) = nanmean(boutdurs);
    statsperflycurr.Z(fly) = nnz(doanalyze);
    statsperflycurr.std(fly) = std(boutdurs,1);
    statsperflycurr.prctiles(:,fly) = prctile(boutdurs,stats_params.prctiles_compute);
    alldata = [alldata,boutdurs]; %#ok<AGROW>
    
  else
    % compute stats
    [statsperflycurr.Z(fly),statsperflycurr.mean(fly),...
      statsperflycurr.std(fly),statsperflycurr.prctiles(:,fly)] = ...
      ComputePerFrameStats(data,doanalyze,...
      'prctiles_compute',stats_params.prctiles_compute);
    alldata = [alldata,data(doanalyze)]; %#ok<AGROW>
    
  end
  
  % we will weight this trajectory by the fraction of frames analyzed
  statsperflycurr.fracframesanalyzed(fly) = nnz(doanalyze_fly) / n;
  statsperflycurr.Zfly(fly) = nnz(doanalyze_fly);

end

if ~isempty(haxcurr),
  for i = 1:numel(hcurr),
    if ~ishandle(hcurr(i)) || hcurr(i) == 0, continue; end
    haxcurr1 = get(hcurr(i),'Parent');
    ylim = get(haxcurr1,'YLim');
    ydata = get(hcurr(i),'YData');
    n3 = numel(ydata)/5;
    set(hcurr(i),'YData',ylim([zeros(n3,2),ones(n3,2),zeros(n3,1)]+1)');
  end
end

function statsperflycurr = ComputeNormStat(statsperflycurr,norm_statsperflycurr)

statsperflycurr.raw_mean = statsperflycurr.mean;
statsperflycurr.norm_mean = norm_statsperflycurr.mean;
statsperflycurr.mean = statsperflycurr.mean - norm_statsperflycurr.mean;
statsperflycurr.raw_Z = statsperflycurr.Z;
statsperflycurr.norm_Z = norm_statsperflycurr.Z;
statsperflycurr.Z = 1./(1./statsperflycurr.raw_Z + 1./statsperflycurr.norm_Z);
statsperflycurr.raw_std = statsperflycurr.std;
statsperflycurr.norm_std = norm_statsperflycurr.std;
statsperflycurr.std = sqrt(statsperflycurr.std.^2 + norm_statsperflycurr.std.^2);
statsperflycurr.raw_prctiles = statsperflycurr.prctiles;
statsperflycurr.norm_prctiles = norm_statsperflycurr.prctiles;
statsperflycurr.prctiles = statsperflycurr.prctiles - norm_statsperflycurr.prctiles;
statsperflycurr.norm_fracframesanalyzed = norm_statsperflycurr.fracframesanalyzed;
statsperflycurr.raw_fracframesanalyzed = statsperflycurr.fracframesanalyzed;
statsperflycurr.fracframesanalyzed = 1./(1./statsperflycurr.raw_fracframesanalyzed+1./statsperflycurr.norm_fracframesanalyzed);
statsperflycurr.raw_Zfly = statsperflycurr.Zfly;
statsperflycurr.norm_Zfly = norm_statsperflycurr.Zfly;
statsperflycurr.Zfly = 1./(1./statsperflycurr.Zfly+1./statsperflycurr.norm_Zfly);

function [data,doanalyze,doanalyze_fly,n,h] = ...
  ComputeDataFly(trx,fly,field,flyconditionname,frameconditionname,flyconditionparams,frameconditionparams,special_cases,varargin)

data = [];

% current field data
if ~ismember(field,special_cases),
  data = trx(fly).(field);
  n = numel(data);
else
  n = trx(fly).nframes;
end

% check that the fly matches the conditions in flyconditionparams
if strcmpi(flyconditionname,'any') && isempty(flyconditionparams),
  doanalyze_fly = true(1,n);
else
  doanalyze_fly = FrameConditionCheck(trx,fly,n,flyconditionparams);
end

% choose frames that match the conditions in frameconditionparams
if strcmpi(frameconditionname,'any'),
  doanalyze_frame = true(1,n);
  h = [];
else
  [doanalyze_frame,h] = FrameConditionCheck(trx,fly,n,frameconditionparams,varargin{:});
end

doanalyze = doanalyze_fly & doanalyze_frame;
