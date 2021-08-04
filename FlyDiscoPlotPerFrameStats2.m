function FlyDiscoPlotPerFrameStats2(expdir,varargin)

version = '0.1';

[analysis_protocol,settingsdir,datalocparamsfilestr,visible,DEBUG,...
  plotstats,plotstim,makestimvideos,plothist,...
  plotflies,plotstimtrajs] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'visible','',...
  'debug',false,...
  'plotstats',true,...
  'plotstim',[],...
  'makestimvideos',[],...
  'plothist',true,...
  'plotflies',true,...
  'plotstimtrajs',true);

if ischar(plotstats),
  plotstats = str2double(plotstats);
  if isnan(plotstats),
    plotstats= true;
  else
    plotstats= plotstats ~= 0;
  end
end
if isempty(visible),
  if DEBUG,
    visible = 'on';
  else
    visible = 'off';
  end
end

%% data locations

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

%% log file

if isfield(dataloc_params,'plotperframestats_logfilestr') && ~DEBUG,
  logfile = fullfile(expdir,dataloc_params.plotperframestats_logfilestr);
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

fprintf(logfid,'\n\n***\nRunning FlyBowlPlotPerFrameStats2 version %s analysis_protocol %s (linked to %s) at %s\n',version,analysis_protocol,real_analysis_protocol,timestamp);

%% load experiment data

statsmatsavename = fullfile(expdir,dataloc_params.statsperframematfilestr);
load(statsmatsavename,'statsperfly','statsperexp');
if plothist,
  histmatsavename = fullfile(expdir,dataloc_params.histperframematfilestr);
  load(histmatsavename,'histperfly','histperexp');
end

%% create the plot directory if it does not exist
figdir = fullfile(expdir,dataloc_params.figdir);
if ~DEBUG && ~exist(figdir,'file'),
  [status,msg] = mkdir(figdir);
  if ~status,
    error('Could not create the figure directory %s:\n%s',figdir,msg);
  end
end

%% load stats params

%statsperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statsperframefeaturesfilestr);
%stats_perframefeatures = ReadStatsPerFrameFeatures2(statsperframefeaturesfile);

statsplotparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statsplotparamsfilestr);
stats_plotparams = ReadStatsPlotParams(statsplotparamsfile);
stimulusplotparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.stimulusplotparamsfilestr);
stimulus_plotparams = ReadParams(stimulusplotparamsfile);

%% load hist params

histperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframefeaturesfilestr);
hist_perframefeatures = ReadStatsPerFrameFeatures2(histperframefeaturesfile);
histperframebinsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframebinsfilestr);
load(histperframebinsfile,'bins');
% statframeconditionsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statframeconditionfilestr);
% frameconditiondict = ReadParams(statframeconditionsfile);

%% read plotting parameters

histplotparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histplotparamsfilestr);
hist_plot_params = ReadHistPlotParams(histplotparamsfile);
conditionplotparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.conditionplotparamsfilestr);
condition_plot_params = ReadConditionParams(conditionplotparamsfile);
[~,expname] = fileparts(expdir);

%% plot means, stds

hfig = [];
hfigperfly = [];

[~,basename] = fileparts(expdir);

statfiles = {};
statperflyfiles = {};

if plotstats,
  groups = fieldnames(stats_plotparams);
  statfiles = cell(size(groups));
  statperflyfiles = cell(size(groups));
    
  for i = 1:numel(groups),
    group = groups{i};
    name = sprintf('%s, %s',basename,group);
    %stathandles = PlotPerFrameStats(stats_plotparams.(group),statsperfly,statsperexp,controlstats,name,'visible',visible,'hfig',i);
    %drawnow;
    % plot a second time to get this to work on the cluster
    stathandles = PlotPerFrameStatsGroup(stats_plotparams.(group),statsperfly,statsperexp,name,'visible',visible,'hfig',hfig,'hfigflies',hfigperfly,'plotparams',condition_plot_params,...
      'plotflies',plotflies);
    if ~DEBUG,
      savename = sprintf('stats_%s.png',group);
      relsavename = fullfile(dataloc_params.figdir,savename);
      savename = fullfile(figdir,savename);
      if exist(savename,'file'),
        delete(savename);
      end
      set(stathandles.hfig,'Units','pixels','Position',stathandles.position);
      save2png(savename,stathandles.hfig);
      statfiles{i} = relsavename;
      hfig = stathandles.hfig;
      
      if plotflies,
        savename = sprintf('stats_perfly_%s.png',group);
        relsavename = fullfile(dataloc_params.figdir,savename);
        savename = fullfile(figdir,savename);
        if exist(savename,'file'),
          delete(savename);
        end
        set(stathandles.hfigflies,'Units','pixels','Position',stathandles.position);
        save2png(savename,stathandles.hfigflies);
        statperflyfiles{i} = relsavename;
        hfigperfly = stathandles.hfigflies;
      end
      
    end
  end
  
end

%% plot stimulus intervals

if isempty(plotstim),
  if isfield(dataloc_params,'indicatordatafilestr'),
    indicatordatafile = fullfile(expdir,dataloc_params.indicatordatafilestr);
    plotstim = exist(indicatordatafile,'file');
  else
    plotstim = false;
  end
end
if isempty(makestimvideos),
  makestimvideos = plotstim;
end

stimfiles = {};
if plotstim || makestimvideos,
  trx = FBATrx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir,...
    'datalocparamsfilestr',datalocparamsfilestr,...
    'maxdatacached',2^30);
  
  %fprintf('Loading trajectories for %s...\n',expdir);
  
  trx.AddExpDir(expdir,'dooverwrite',false,'openmovie',false);
  ind = trx.getIndicatorLED(1);
  non = numel(ind.starton);
  if non == 0,
    plotstim = false;
    makestimvideos = false;
  end
end

hfigperfly = [];

if plotstim,
  stimfiles = cell(non,numel(stimulus_plotparams.features));
  stimperflyfiles = cell(non,numel(stimulus_plotparams.features));
  
  for fi = 1:numel(stimulus_plotparams.features),
    field = stimulus_plotparams.features{fi};
    if isfield(bins,field),
      ylim = bins.(field).edges_linear([1,end]);
    else
      data = trx.(field);
      data = cat(2,data{:});
      ylim = prctile(data(~isnan(data)),[0,100]+[1,-1]*stimulus_plotparams.ylim_prctile);
    end
    for ion = 1:non,
      [hfigcurr,hfigperflycurr] = PlotStimulusInterval(trx,field,ion,basename,...
        'visible',visible,'prestim',stimulus_plotparams.prestim,...
        'poststim',stimulus_plotparams.poststim,...
        'ylim',ylim,'plotflies',plotflies,...
        'hfig',hfig,...
        'hfigfly',hfigperfly,...
        'maxnflies',stimulus_plotparams.maxnflies);
      
      if ~DEBUG,
        
        if ishandle(hfigcurr),
          savename = sprintf('stimulus_%d_%s.png',ion,field);
          relsavename = fullfile(dataloc_params.figdir,savename);
          savename = fullfile(figdir,savename);
          if exist(savename,'file'),
            delete(savename);
          end
          
          set(hfigcurr,'Units','pixels');
          save2png(savename,hfigcurr);
          stimfiles{ion,fi} = relsavename;
          hfig = hfigcurr;
        else
          warning('Could not save stimulus plot for period %d, field %s, figure handle invalid',ion,field);
        end
        
        if plotflies,
          if ishandle(hfigperflycurr),
            savename = sprintf('stimulusperfly_%d_%s.png',ion,field);
            relsavename = fullfile(dataloc_params.figdir,savename);
            savename = fullfile(figdir,savename);
            if exist(savename,'file'),
              delete(savename);
            end
            
            set(hfigperflycurr,'Units','pixels');
            save2png(savename,hfigperflycurr);
            stimperflyfiles{ion,fi} = relsavename;
            hfigperfly = hfigperflycurr;
          else
            warning('Could not save per-fly stimulus plot for period %d, field %s, figure handle invalid',ion,field);
          end
        end
        
      end
      
    end
  end
end

%% make plots of trajectories at start of stimuli for individual flies

if plotstimtrajs,
  [hfigcurr] = PlotStimulusOnsetTrajs(trx,...
    'prestim',stimulus_plotparams.traj_prestim,'poststim',stimulus_plotparams.traj_poststim,...
    'minboxwidth',stimulus_plotparams.traj_minboxwidth,'boxborder',stimulus_plotparams.traj_boxborder,...
    'downsample',stimulus_plotparams.traj_downsample,...
    'nfliesplot',stimulus_plotparams.traj_nfliesplot,...
    'nperiodsplot',stimulus_plotparams.traj_nperiodsplot,...
    'maxnflies',stimulus_plotparams.maxnflies,...
    'hfig',hfig,...
    'visible',visible);
  if ~DEBUG,
    
    if ishandle(hfigcurr),
      savename = 'stimulus_traj.png';
      relsavename = fullfile(dataloc_params.figdir,savename);
      savename = fullfile(figdir,savename);
      if exist(savename,'file'),
        delete(savename);
      end
      
      set(hfigcurr,'Units','pixels');
      save2png(savename,hfigcurr);
      stimtrajfile = relsavename;
      hfig = hfigcurr;
    else
      warning('Could not save stimulus trajectory plot');
    end
    
  end
else
  stimtrajfile = '';
end
%% make videos at start of stimuli for individual flies

videofiles = {};
ions = [];
fliesplot = [];

if makestimvideos,
  
  stimulusvideoparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.stimulusvideoparamsfilestr);
  videoparams = ReadParams(stimulusvideoparamsfile);
  videobasefile = fullfile(expdir,dataloc_params.figdir,dataloc_params.stimulusvideofilestr);
  moviefile = fullfile(expdir,dataloc_params.moviefilestr);
  if videoparams.nperiodsplot < 0,
    videoparams.nperiodsplot = [];
  end
  if videoparams.nfliesplot < 0,
    videoparams.nfliesplot = [];
  end
  
  [videofiles,ions,fliesplot] = MakeStimulusOnsetMovies(moviefile,trx,videobasefile,...
    'prestim',videoparams.prestim,'poststim',videoparams.poststim,...
    'minboxwidth',videoparams.minboxwidth,'boxborder',videoparams.boxborder,...
    'downsample',videoparams.downsample,...
    'nfliesplot',videoparams.nfliesplot,...
    'nperiodsplot',videoparams.nperiodsplot,...
    'visible',visible);
  for i = 1:numel(videofiles),
    videofiles{i} = videofiles{i}(numel(expdir)+2:end);
  end
  
end

%% plot histograms

histfiles = {};
if plothist,

  ngroups = numel(hist_plot_params.groups);
  histfiles = cell(1,ngroups);
  
  for i = 1:ngroups,
  
    group = hist_plot_params.groups(i);
    [handles,didplot] = PlotPerFrameHistsGroup(group,hist_perframefeatures,...
      histperexp,histperfly,...
      bins,hist_plot_params,condition_plot_params,expname,...
      'visible',visible,...
      'hfig',hfig);
  
    if ~didplot,
      continue;
    end

    if ~DEBUG,
      savename = sprintf('hist_%s.png',group.name);
      relsavename = fullfile(dataloc_params.figdir,savename);
      savename = fullfile(figdir,savename);      
      if exist(savename,'file'),
        delete(savename);
      end
      set(handles.hfig,'Units','pixels','Position',handles.position);
      save2png(savename,handles.hfig);
      histfiles{i} = relsavename;
    end
  end

end

%% output summary webpage with all plots embedded
if ~DEBUG,
  statswebpagefile = fullfile(expdir,dataloc_params.statswebpagefilestr);
  if exist(statswebpagefile,'file'),
    delete(statswebpagefile);
  end
  MakePerFrameStatsWebpage(statswebpagefile,basename,groups,statfiles,...
    statperflyfiles,stimulus_plotparams.features,stimfiles,stimperflyfiles,...
    ions,fliesplot,stimtrajfile,videofiles,...
    histfiles,{hist_plot_params.groups.name});
end

%% save info to a mat file

filename = fullfile(expdir,dataloc_params.plotperframestatsinfomatfilestr);
fprintf(logfid,'Saving debug info to file %s...\n',filename);
ppfsinfo = struct;
ppfsinfo.params = struct;
ppfsinfo.params.stats_plotparams = stats_plotparams;
ppfsinfo.params.stimulus_plotparams = stimulus_plotparams;
ppfsinfo.params.hist_plot_params = hist_plot_params;
ppfsinfo.params.condition_plot_params = condition_plot_params;
ppfsinfo.params.videoparams = videoparams;

ppfsinfo.statsmatsavename = statsmatsavename;
if plothist,
  ppfsinfo.histmatsavename = histmatsavename;
end

ppfsinfo.version = version;
ppfsinfo.analysis_protocol = analysis_protocol;
ppfsinfo.linked_analysis_protocol = real_analysis_protocol;
ppfsinfo.timestamp = timestamp;

if exist(filename,'file'),
  try %#ok<TRYNC>
    delete(filename);
  end
end
if ~DEBUG,
  try
    save(filename,'-struct','ppfsinfo');
  catch ME,
    warning('Could not save information to file %s: %s',filename,getReport(ME));
  end
end

if ~DEBUG,
  close all;
end

