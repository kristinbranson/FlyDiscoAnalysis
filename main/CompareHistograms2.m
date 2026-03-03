function CompareHistograms2(expdirs,outdir,varargin)

[analysis_protocol,settingsdir,datalocparamsfilestr,visible,controldatadirstr] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'visible','off',...
  'controldatadirstr','current');

nexps = numel(expdirs);

%% data locations

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

%% load experiment data

statsdata = cell(1,nexps);
histdata = cell(1,nexps);
for expi = 1:nexps,
  expdir = expdirs{expi};
  statsmatsavename = fullfile(expdir,dataloc_params.statsperframematfilestr);
  statsdata{expi} = load(statsmatsavename);
  histmatsavename = fullfile(expdir,dataloc_params.histperframematfilestr);
  histdata{expi} = load(histmatsavename);
end

%% create the plot directory if it does not exist
figdir = fullfile(outdir,dataloc_params.figdir);
if ~exist(figdir,'file'),
  [status,msg] = mkdir(figdir);
  if ~status,
    error('Could not create the figure directory %s:\n%s',figdir,msg);
  end
end

%% load stats params

statsperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statsperframefeaturesfilestr);
stats_perframefeatures = ReadStatsPerFrameFeatures(statsperframefeaturesfile);

%% load hist params

histperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframefeaturesfilestr);
hist_perframefeatures = ReadHistPerFrameFeatures(histperframefeaturesfile);
histperframebinsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframebinsfilestr);
load(histperframebinsfile,'bins');

%% read plotting parameters

histplotparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histplotparamsfilestr);
hist_plot_params = ReadParams(histplotparamsfile);
expnames = cell(1,nexps);
for expi = 1:nexps,
  [~,expnames{expi}] = fileparts(expdirs{expi});
end

%% get control data

% if ~isempty(controldatadirstr),
%   
%   controldatadir = fullfile(dataloc_params.pBDPGAL4Ustatsdir,controldatadirstr);
%   controlstatsname = fullfile(controldatadir,dataloc_params.statsperframematfilestr);
%   controlstats = load(controlstatsname);
%   controlhistname = fullfile(controldatadir,dataloc_params.histperframematfilestr);
%   controlhist = load(controlhistname);
%   
% else
%   
%   controlstats = [];
%   controlhist = [];
%   
% end

%% plot means, stds

% for expi = 1:nexps,
%   expdir = expdirs{expi};
%   [~,basename] = fileparts(expdir);
%   stathandles = PlotPerFrameStats(stats_perframefeatures,statsperfly,statsperexp,controlstats,basename,'visible',visible);
%   drawnow;
% end
% savename = sprintf('stats.png');
% savename = fullfile(figdir,savename);
% if exist(savename,'file'),
%   delete(savename);
% end
% save2png(savename,stathandles.hfig);

%% plot histograms

hist_fields = unique({hist_perframefeatures.field});
linestyles = {'-',':','--','-.'};
for i = 1:numel(hist_fields),
  
  field = hist_fields{i};
  hax = [];
  handles = cell(1,nexps);
  
  for expi = nexps:-1:1,
    
    linestyle = linestyles{mod(expi-1,numel(linestyles))+1};
    
    if isfield(histdata{expi},'meanhistperexp'),
      handles{expi} = PlotPerExpHists(field,hist_perframefeatures,...
        histdata{expi}.meanhistperexp,histdata{expi}.histperexp,...
        bins.(field),hist_plot_params,'',...
        'visible',visible,'linestyle',linestyle,'stdstyle','errorbar','hax',hax);
    else
      handles{expi} = PlotPerFrameHists(field,hist_perframefeatures,...
        histdata{expi}.histperexp,histdata{expi}.histperfly,...
        bins.(field),hist_plot_params,'',...
        'visible',visible,'linestyle',linestyle,'stdstyle','errorbar','hax',hax);
    end
    hax = handles{expi}.hax;
  end
  
  % fix legend
  hs = nan(1,nexps);
  for expi = 1:nexps,
    hs(expi) = handles{expi}.htype(1);
  end
  s = get(handles{1}.hleg,'String');
  legend([hs,handles{1}.htype],[expnames,s],'Parent',handles{1}.hfig,'Interpreter','none','Location','best');
  
  drawnow;
  savename = sprintf('hist_%s.png',hist_fields{i});
  savename = fullfile(figdir,savename);
  if exist(savename,'file'),
    delete(savename);
  end
  save2png(savename,handles{end}.hfig);
  
end

close all;