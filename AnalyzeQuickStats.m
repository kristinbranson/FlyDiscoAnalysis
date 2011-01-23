function quickstats_stats = AnalyzeQuickStats(varargin)

%% parameters

[experiment_params,UFMFDiagnosticsFileStr,...
  ufmfstream_nbins,savename] = myparse(varargin,...
  'experiment_params',{},...
  'UFMFDiagnosticsFileStr','ufmf_diagnostics.txt',...
  'ufmfstream_nbins',100,...
  'savename',['QuickStats_Stats_',datestr(now,30),'.mat']);

%% experiment directories


[expdirs,expdir_reads,expdir_writes,experiments] = ...
  getExperimentDirs(experiment_params{:});

nexpdirs = numel(expdirs);

%% compute stats of quick stats and ufmf diagnostics

meanquickstats = struct;
stdquickstats = struct;
meanufmfstats = struct;
stdufmfstats = struct;
%streamcounts = struct;
%streamlengths = struct;
streamstats = struct;

for i = 1:nexpdirs,
  
  % load the quick stats
  quickstats = loadQuickStats(expdir_reads{i});
  
  % read ufmf diagnostics
  UFMFDiagnosticsFileName = fullfile(expdir_reads{i},UFMFDiagnosticsFileStr);
  [UFMFStats,success1,errmsg] = readUFMFDiagnostics(UFMFDiagnosticsFileName);
  if ~success1,
    error('Could not read ufmf diagnostics file %s: %s',UFMFDiagnosticsFileName,errmsg);
  end
  
  % timestamp offset
  UFMFStats.stream.timestamp = UFMFStats.stream.timestamp - UFMFStats.stream.timestamp(1);

  % first experiment, don't know field names
  if i == 1,
    % initialize means, stds for quickstats
    meanquickstats = quickstats;
    fns = fieldnames(quickstats);
    for j = 1:numel(fns),
      stdquickstats.(fns{j}) = quickstats.(fns{j}).^2;
    end
    
    % fieldnames for stream, summary
    streamfns = fieldnames(UFMFStats.stream);
    summaryfns = fieldnames(UFMFStats.summary);
    
    % store all stream stats cuz we don't know where to bin yet
    streamstats = UFMFStats.stream;
    
    % initialize means, stds for summary
    meanufmfstats.summary = UFMFStats.summary;
    for j = 1:numel(summaryfns),
      stdufmfstats.summary.(summaryfns{j}) = UFMFStats.summary.(summaryfns{j}).^2;
    end
    
  % not the first experiment
  else
    
    % increment mean, std quick stats
    for j = 1:numel(fns),
      meanquickstats.(fns{j}) = meanquickstats.(fns{j}) + quickstats.(fns{j});
      stdquickstats.(fns{j}) = stdquickstats.(fns{j}) + quickstats.(fns{j}).^2;
    end
    
    % store stream stats
    streamstats = structarrayset(streamstats,i,UFMFStats.stream);
    
    % increment means, stds for summary
    for j = 1:numel(summaryfns),
      meanufmfstats.summary.(summaryfns{j}) = meanufmfstats.summary.(summaryfns{j}) + ...
        UFMFStats.summary.(summaryfns{j});
      stdufmfstats.summary.(summaryfns{j}) = stdufmfstats.summary.(summaryfns{j}) + ...
        UFMFStats.summary.(summaryfns{j}).^2;
    end
    
  end
end

% normalize quick stats
for j = 1:numel(fns),
  meanquickstats.(fns{j}) = meanquickstats.(fns{j}) ./ nexpdirs;
  stdquickstats.(fns{j}) = sqrt(stdquickstats.(fns{j}) ./ nexpdirs - meanquickstats.(fns{j}).^2);
end

% normalize summary stats
for j = 1:numel(summaryfns),
  meanufmfstats.summary.(summaryfns{j}) = meanufmfstats.summary.(summaryfns{j}) ./ nexpdirs;
  stdufmfstats.summary.(summaryfns{j}) = sqrt(stdufmfstats.summary.(summaryfns{j}) ./ nexpdirs - meanufmfstats.summary.(summaryfns{j}).^2);
end

% compute stats for streams
dt = [streamstats.timestamp];
mindt = min(dt);
maxdt = max(dt);
dt_edges = linspace(mindt,maxdt,ufmfstream_nbins+1);
[~,idx] = histc(dt,dt_edges);
for i = 1:numel(streamfns),
  fn = streamfns{i};
  tmp = [streamstats.(fn)];
  meanufmfstats.stream.(fn) = nan(1,ufmfstream_nbins);
  stdufmfstats.stream.(fn) = nan(1,ufmfstream_nbins);
  for j = 1:ufmfstream_nbins,
    meanufmfstats.stream.(fn)(j) = nanmean(tmp(idx==j));
    stdufmfstats.stream.(fn)(j) = nanstd(tmp(idx==j));
  end
end

%% Store in a format easy to load into computeQuickStats

quickstats_stats = struct;
quickstats_stats.IntensityHistMu = meanquickstats.BkgdIntensityHist_frac;
quickstats_stats.IntensityHistSig = stdquickstats.BkgdIntensityHist_frac;
quickstats_stats.ScanLineMu = [];
quickstats_stats.ScanLineSig = [];
i = 1;
while true,
  fn = sprintf('BkgdScanLine_intensities_%d',i);
  if ~isfield(meanquickstats,fn),
    break;
  end
  quickstats_stats.ScanLineMu(i,:) = meanquickstats.(fn);
  quickstats_stats.ScanLineSig(i,:) = stdquickstats.(fn);
  i = i + 1;
end
quickstats_stats.UFMFStreamMu = struct;
quickstats_stats.UFMFStreamSig = struct;
fns = fieldnames(meanufmfstats.stream);
for i = 1:numel(fns),
  quickstats_stats.UFMFStreamMu.(fns{i}) = meanufmfstats.stream.(fns{i});
  quickstats_stats.UFMFStreamSig.(fns{i}) = stdufmfstats.stream.(fns{i});
end
quickstats_stats.UFMFSummaryMu = struct;
quickstats_stats.UFMFSummarySig = struct;
fns = fieldnames(meanufmfstats.summary);
for i = 1:numel(fns),
  quickstats_stats.UFMFSummaryMu.(fns{i}) = meanufmfstats.summary.(fns{i});
  quickstats_stats.UFMFSummarySig.(fns{i}) = stdufmfstats.summary.(fns{i});
end

save(savename,'-struct','quickstats_stats');