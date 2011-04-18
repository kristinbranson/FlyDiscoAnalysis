function FlyBowlPlotPerFrameStats(expdir,varargin)

[analysis_protocol,settingsdir,datalocparamsfilestr,visible,controldatadirstr] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'visible','off',...
  'controldatadirstr','current');

%% data locations

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

%% load experiment data

statsmatsavename = fullfile(expdir,dataloc_params.statsperframematfilestr);
load(statsmatsavename,'statsperfly','statsperexp');
histmatsavename = fullfile(expdir,dataloc_params.histperframematfilestr);
load(histmatsavename,'histperfly','histperexp');

%% create the plot directory if it does not exist
figdir = fullfile(expdir,dataloc_params.figdir);
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
[tmp,expname] = fileparts(expdir);

%% get control data

if ~isempty(controldatadirstr),
  
  controldatadir = fullfile(dataloc_params.pBDPGAL4Ustatsdir,controldatadirstr);
  controlstatsname = fullfile(controldatadir,dataloc_params.statsperframematfilestr);
  controlstats = load(controlstatsname);
  controlhistname = fullfile(controldatadir,dataloc_params.histperframematfilestr);
  controlhist = load(controlhistname);
  
else
  
  controlstats = [];
  controlhist = [];
  
end

%% plot means, stds

[tmp,basename] = fileparts(expdir);
stathandles = PlotPerFrameStats(stats_perframefeatures,statsperfly,statsperexp,controlstats,basename,'visible',visible);
drawnow;  
savename = sprintf('stats.png');
savename = fullfile(figdir,savename);
if exist(savename,'file'),
  delete(savename);
end
save2png(savename,stathandles.hfig);

%% plot histograms

hist_fields = unique({hist_perframefeatures.field});
for i = 1:numel(hist_fields),
  
  field = hist_fields{i};

  if ~isempty(controlhist),
    if isfield(controlhist,'meanhistperexp'),
      handles_control = PlotPerExpHists(field,hist_perframefeatures,...
        controlhist.meanhistperexp,controlhist.histperexp,...
        bins.(field),hist_plot_params,expname,...
        'visible',visible,'linestyle',':','stdstyle','errorbar');
    else
      handles_control = PlotPerFrameHists(field,hist_perframefeatures,...
        controlhist.meanhistperfly,controlhist.histperfly,...
        bins.(field),hist_plot_params,expname,...
        'visible',visible,'linestyle',':','stdstyle','errorbar');
    end
    hax = handles_control.hax;
  else
    hax = [];
  end
  
  handles = PlotPerFrameHists(field,hist_perframefeatures,...
    histperexp,histperfly,...
    bins.(field),hist_plot_params,expname,...
    'visible',visible,...
    'hax',hax);
  
  if ~isempty(controlhist),
    % fix legend
    s = get(handles.hleg,'String');
    legend([handles_control.htype(1),handles.htype],[{'control'},s],'Parent',handles.hfig,'Interpreter','none');
  end
  
  drawnow;
  savename = sprintf('hist_%s.png',hist_fields{i});
  savename = fullfile(figdir,savename);
  if exist(savename,'file'),
    delete(savename);
  end
  save2png(savename,handles.hfig);
  
end

close all;