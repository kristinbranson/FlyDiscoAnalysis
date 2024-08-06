% FlyDiscoPlotPerFrameStats(expdir,varargin)
% Optional arguments:
% forcecompute: Whether the whole stage has been forced (usually used in the
%   context of FlyDiscoPipeline())
%     false -> The stage has not been forced, so leave pre-existing outputs in
%              place
%     true -> The stage *has* been forced, so delete any pre-existing HTML
%             output file and the output folder before running anything
% Default value: false
% plottimeseries: Whether to plot mean stats. 
%     0 -> Do not plot
%     1 -> Plot if files do not exist
%     2 -> Replot whether or not files exist
% Default value 2. 
% plotstats: Whether to plot mean stats. 
%     0 -> Do not plot
%     1 -> Plot if files do not exist
%     2 -> Replot whether or not files exist
% Default value 2. 
% plotstim: Whether to plot stimulus onset time series intervals
%     0 -> Do not plot
%     1 -> Plot if files do not exist
%     2 -> Replot whether or not files exist
%     [] -> 0 if no indicator file, 2 if indicator file
% Default value [].
% plotstimtrajs: Whether to plot stimulus onset trajectories
%     0 -> Do not plot
%     1 -> Plot if files do not exist
%     2 -> Replot whether or not files exist
%     [] -> 0 if no indicator file, 2 if indicator file
% Default value [].
% plothist: Whether to plot histograms of statistics
%     0 -> Do not plot
%     1 -> Plot if files do not exist
%     2 -> Replot whether or not files exist
% Default value 2.
% makestimvideos: Whether to make gifs of stimulus onset periods. 
%     0 -> Do not plot
%     1 -> Plot if files don't exist
%     2 -> Replot whether or not files exist. 
%     [] -> Inherit value from plotstimtrajs
% Default value [].
% plotflies: Whether to plot per-fly mean statistics. Only relevant if
% plotstats > 0. Binary false/true value. 
% Default value true.
% analysis_protocol: Default = 'current'
% settingsdir: Default = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings'
% datalocparamsfilestr: Default = 'dataloc_params.txt'
% visible: Whether to show the figures while plotting. If empty, then this
% is true when debug parameter is true and false otherwise. Default value
% ''.
% debug: Whether to just plot everything (true) or to save plots to
% files. 
% verbose: 0-2 value indicating how much to print to stdout. Default = 1.

function FlyDiscoPlotPerFrameStats(expdir,varargin)

version = '0.1';

[analysis_protocol,settingsdir,datalocparamsfilestr,visible,DEBUG,...
  forcecompute,plottimeseries,plotstats,plotstim,makestimvideos,plothist,...
  plotflies,plotstimtrajs,verbose] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'visible','',...
  'debug',false,...
  'forcecompute', false, ...
  'plottimeseries',2,...
  'plotstats',2,...
  'plotstim',[],...
  'makestimvideos',[],...
  'plothist',2,...
  'plotflies',true,...
  'plotstimtrajs',[],...
  'verbose',1);

if ischar(plotstats),
  plotstats = str2double(plotstats);
  if isnan(plotstats),
    plotstats= 2;
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

% Paths to outputs
figdir = fullfile(expdir, dataloc_params.figdir) ;
statswebpagefile = fullfile(expdir,dataloc_params.statswebpagefilestr);

% If we're forcing everything, delete the outputs before we do anything else
if forcecompute ,
    ensure_file_or_folder_does_not_exist(figdir) ;
    ensure_file_does_not_exist(statswebpagefile) ;
end

% log file is always stdout
% if isfield(dataloc_params,'plotperframestats_logfilestr') && ~DEBUG,
%   logfile = fullfile(expdir,dataloc_params.plotperframestats_logfilestr);
%   logfid = fopen(logfile,'a');
%   if logfid < 1,
%     warning('Could not open log file %s\n',logfile);
%     logfid = 1;
%   end
% else
logfid = 1;
% end

timestamp = datestr(now,'yyyymmddTHHMMSS');
real_analysis_protocol = GetRealAnalysisProtocol(analysis_protocol,settingsdir);

%fprintf(logfid,'\n\n***\nRunning FlyDiscoPlotPerFrameStats2 version %s analysis_protocol %s (linked to %s) at %s\n',version,analysis_protocol,real_analysis_protocol,timestamp);

%% load experiment data

statsmatsavename = fullfile(expdir,dataloc_params.statsperframematfilestr);
load(statsmatsavename,'statsperfly','statsperexp');
if plothist,
  histmatsavename = fullfile(expdir,dataloc_params.histperframematfilestr);
  load(histmatsavename,'histperfly','histperexp');
end

%% create the plot directory if it does not exist
if ~DEBUG && ~exist(figdir,'file'),
  [status,msg] = mkdir(figdir);
  if ~status,
    error('Could not create the figure directory %s:\n%s',figdir,msg);
  end
end

%% load stats params

%statsperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statsperframefeaturesfilestr);
%stats_perframefeatures = ReadStatsPerFrameFeatures2(statsperframefeaturesfile);
timeseriesplotparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.timeseriesplotparamsfilestr);
timeseries_params = ReadParams(timeseriesplotparamsfile);
statsplotparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statsplotparamsfilestr);
stats_plotparams = ReadStatsPlotParams(statsplotparamsfile);
stimulusplotparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.stimulusplotparamsfilestr);
stimulus_plotparams = ReadStimulusPlotParams(stimulusplotparamsfile);  

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

%% plot timeseries
% add checks for whether this run
hfig = [];
% plottimeseries = true;
timeseriesfiles = {};
if plottimeseries
    if verbose > 0
        fprintf('Plotting time series...\n')
    end

%     timeseries_params = ReadParams('/groups/branson/home/robiea/Code_versioned/FlyDiscoAnalysis/settings-internal/20230907_flybubble_LED_VNC2_hack/timeseries_perframefeatures.txt');
    trx = FBATrx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir,...
        'datalocparamsfilestr',datalocparamsfilestr,...
        'maxdatacached',2^30);
    trx.AddExpDir(expdir,'dooverwrite',false,'openmovie',false);
    % load indicator data
    ind = trx.getIndicatorLED(1);
    non = numel(ind.starton);
    % check if there are stim periods in indicator data 
    % TODO use this check
    if non == 0,
        plottimeseries = 0;
    end 

    % fix for just 1 feature
    if ischar(timeseries_params.features)
        timeseries_params.features = {timeseries_params.features};
    end
    timeseriesfiles = cell(size(timeseries_params.features));
    for fns = 1:numel(timeseries_params.features)    
        field = timeseries_params.features{fns};
        savename = sprintf('timeseries_%s.png',field);
        relsavename = fullfile(dataloc_params.figdir,savename);
        savename = fullfile(figdir,savename);
        if ~DEBUG &&  plottimeseries <2 && exist(savename,'file')
            timeseriesfiles{fns} = relsavename;
            continue;
        end
        if verbose > 1,
            fprintf('Plotting time series for %s\n',fns);
        end

        inputparams = {};
        if isfield(timeseries_params,'bin')
            inputparams = {'bin', timeseries_params.bin};
        end
        if isfield(timeseries_params,'perctimebin')
            inputparams = [inputparams,'percbin',timeseries_params.perctimebin];
        end
        tshandle = PlotPerFrameTimeSeries(trx,ind,field,'visible',visible,'hfig',hfig,inputparams{:});

        if ~DEBUG,
            if exist(savename,'file'),
                delete(savename);
            end
            set(tshandle,'Units','pixels','Position',tshandle.Position)
            save2png(savename,tshandle);
            timeseriesfiles{fns} = relsavename;
            hfig = tshandle;

        end
    end
end
%% plot means, stds

hfig = [];
hfigperfly = [];

[~,basename] = fileparts(expdir);

statfiles = {};
statperflyfiles = {};

if plotstats,
  
  if verbose > 0,
    fprintf('Plotting mean stats...\n');
  end
  
  groups = fieldnames(stats_plotparams);
  statfiles = cell(size(groups));
  statperflyfiles = cell(size(groups));
    
  for i = 1:numel(groups),
    group = groups{i};
    %name = sprintf('%s, %s',basename,group);
    
    savename = sprintf('stats_%s.png',group);
    relsavename = fullfile(dataloc_params.figdir,savename);
    savename = fullfile(figdir,savename);
    savename_fly = sprintf('stats_perfly_%s.png',group);
    relsavename_fly = fullfile(dataloc_params.figdir,savename_fly);
    savename_fly = fullfile(figdir,savename_fly);

    if ~DEBUG && plotstats < 2 && exist(savename,'file') && ...
        (~plotflies || exist(savename_fly,'file')),
      statfiles{i} = relsavename;
      if plotflies,
        statperflyfiles{i} = relsavename_fly;
      end
      continue;
    end
    if verbose > 1,
      fprintf('Plotting per-frame stats for group %s\n',group);
    end
    stathandles = PlotPerFrameStatsGroup(stats_plotparams.(group),statsperfly,statsperexp,group,'visible',visible,'hfig',hfig,'hfigflies',hfigperfly,'plotparams',condition_plot_params,...
      'plotflies',plotflies);
    if ~DEBUG,
      if exist(savename,'file'),
        delete(savename);
      end
      set(stathandles.hfig,'Units','pixels','Position',stathandles.position);
      save2png(savename,stathandles.hfig);
      statfiles{i} = relsavename;
      hfig = stathandles.hfig;
      
      if plotflies,
        if exist(savename_fly,'file'),
          delete(savename_fly);
        end
        set(stathandles.hfigflies,'Units','pixels','Position',stathandles.position);
        save2png(savename_fly,stathandles.hfigflies);
        statperflyfiles{i} = relsavename_fly;
        hfigperfly = stathandles.hfigflies;
      end
      
    end
  end
  
end

%% plot stimulus intervals

if isempty(plotstim),
  if isfield(dataloc_params,'indicatordatafilestr'),
    indicatordatafile = fullfile(expdir,dataloc_params.indicatordatafilestr);
    if exist(indicatordatafile,'file'),
      plotstim = 2;
    else
      plotstim = 0;
    end
  else
    plotstim = 0;
  end
end
if isempty(plotstimtrajs),
  if isfield(dataloc_params,'indicatordatafilestr'),
    indicatordatafile = fullfile(expdir,dataloc_params.indicatordatafilestr);
    if exist(indicatordatafile,'file'),
      plotstimtrajs = 2;
    else
      plotstimtrajs = 0;
    end
  else
    plotstimtrajs = 0;
  end
end
if isempty(makestimvideos),
  makestimvideos = plotstim;
end

stimfiles = {};
if plotstim > 0 || makestimvideos > 0,
  trx = FBATrx('analysis_protocol',analysis_protocol,'settingsdir',settingsdir,...
    'datalocparamsfilestr',datalocparamsfilestr,...
    'maxdatacached',2^30);
  
  %fprintf('Loading trajectories for %s...\n',expdir);
  
  trx.AddExpDir(expdir,'dooverwrite',false,'openmovie',false);
  ind = trx.getIndicatorLED(1);
  non = numel(ind.starton);
  if non == 0,
    plotstim = 0;
    makestimvideos = 0;
  end

  if ~isfield(stimulus_plotparams,'stimsets'),
    stimulus_plotparams.stimsets = cell(non,2);
    for i = 1:non,
      stimulus_plotparams.stimsets{i,1} = num2str(i);
      stimulus_plotparams.stimsets{i,2} = i;
    end
  end
  nsetson = size(stimulus_plotparams.stimsets,1);


end

hfigperfly = [];

if plotstim > 0,
  
  if verbose > 0,
    fprintf('Plotting per-frame features at stimulus onset...\n');
  end

  stimfiles = cell(nsetson,numel(stimulus_plotparams.features));
  stimperflyfiles = cell(nsetson,numel(stimulus_plotparams.features));
  
  ind = trx.getIndicatorLED(1);
  setstarts = cellfun(@(x) min(x), stimulus_plotparams.stimsets(:,2));
  fliesplot = ChooseFliesPlot(trx,ind,setstarts,stimulus_plotparams.maxnflies);
  
  for fi = 1:numel(stimulus_plotparams.features),
    field = stimulus_plotparams.features{fi};
    if isfield(bins,field),
      ylim = bins.(field).edges_linear([1,end]);
    else
      data = trx.(field);
      data = cat(2,data{:});
      ylim = prctile(data(~isnan(data)),[0,100]+[1,-1]*stimulus_plotparams.ylim_prctile);
    end
    for seti = 1:nsetson,
      setname = stimulus_plotparams.stimsets{seti,1};
      ions = stimulus_plotparams.stimsets{seti,2};
      savename = sprintf('stimulus_%s_%s.png',setname,field);
      relsavename = fullfile(dataloc_params.figdir,savename);
      savename = fullfile(figdir,savename);
      savename_fly = sprintf('stimulusperfly_%s_%s.png',setname,field);
      relsavename_fly = fullfile(dataloc_params.figdir,savename_fly);
      savename_fly = fullfile(figdir,savename_fly);
      
      if ~DEBUG && plotstim < 2 && exist(savename,'file') && ...
          (~plotflies || exist(savename_fly,'file')),
        stimfiles{seti,fi} = relsavename;
        if plotflies,
          stimperflyfiles{seti,fi} = relsavename_fly;
        end
        continue;
      end
      
      if verbose > 1,
        fprintf('Plotting per-frame features at stimulus onset, feature %s, interval %d',field,ion);
      end
      [hfigcurr,hfigperflycurr] = PlotStimulusInterval(trx,field,ions,basename,setname,...
        'visible',visible,'prestim',stimulus_plotparams.prestim,...
        'poststim',stimulus_plotparams.poststim,...
        'ylim',ylim,'plotflies',plotflies,...
        'hfig',hfig,...
        'hfigfly',hfigperfly,...
        'fliesplot',fliesplot);
      
      if ~DEBUG,
        
        if ishandle(hfigcurr),
          if exist(savename,'file'),
            delete(savename);
          end
          
          set(hfigcurr,'Units','pixels');
          save2png(savename,hfigcurr);
          stimfiles{seti,fi} = relsavename;
          hfig = hfigcurr;
        else
          warning('Could not save stimulus plot for period %d, field %s, figure handle invalid',ion,field);
        end
        
        if plotflies,
          if ishandle(hfigperflycurr),
            if exist(savename_fly,'file'),
              delete(savename_fly);
            end
            
            set(hfigperflycurr,'Units','pixels');
            save2png(savename_fly,hfigperflycurr);
            stimperflyfiles{seti,fi} = relsavename_fly;
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

if plotstimtrajs > 0,
  
  savename = 'stimulus_traj.png';
  relsavename = fullfile(dataloc_params.figdir,savename);
  savename = fullfile(figdir,savename);
  
  if verbose > 0,
    fprintf('Plotting stimulus onset trajectories...\n');
  end
  
  if ~DEBUG && plotstimtrajs < 2 && exist(savename,'file'),
    stimtrajfile = relsavename;
  else
    [hfigcurr] = PlotStimulusOnsetTrajs(trx,...
      'prestim',stimulus_plotparams.traj_prestim,'poststim',stimulus_plotparams.traj_poststim,...
      'minboxwidth',stimulus_plotparams.traj_minboxwidth,'boxborder',stimulus_plotparams.traj_boxborder,...
      'downsample',stimulus_plotparams.traj_downsample,...
      'nfliesplot',stimulus_plotparams.traj_nfliesplot,...
      'nperiodsplot',stimulus_plotparams.traj_nperiodsplot,...
      'stimsets',stimulus_plotparams.stimsets,...
      'maxnflies',stimulus_plotparams.maxnflies,...
      'hfig',hfig,...
      'visible',visible);
%       'plotstd',stimulus_plotparams.plotstd,...
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
  end
else
  stimtrajfile = '';
end
%% make videos at start of stimuli for individual flies

videofiles = {};
videoions = [];
fliesplot = [];

if makestimvideos,
  
  if verbose > 0,
    fprintf('Making stimulus onset videos...\n');
  end
  
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
  
  [videofiles,videoions,fliesplot] = MakeStimulusOnsetMovies(moviefile,trx,videobasefile,...
    'prestim',videoparams.prestim,'poststim',videoparams.poststim,...
    'minboxwidth',videoparams.minboxwidth,'boxborder',videoparams.boxborder,...
    'downsample',videoparams.downsample,...
    'nfliesplot',videoparams.nfliesplot,...
    'nperiodsplot',videoparams.nperiodsplot,...
    'stimsets',stimulus_plotparams.stimsets,...
    'visible',visible,...
    'force',makestimvideos>1);
  for i = 1:numel(videofiles),
    videofiles{i} = videofiles{i}(numel(expdir)+2:end);
  end
  
end

%% plot histograms

histfiles = {};
if plothist,

  ngroups = numel(hist_plot_params.groups);
  histfiles = cell(1,ngroups);
  if verbose > 0,
    fprintf('Plotting histograms...\n');
  end
  
  for i = 1:ngroups,
  
    group = hist_plot_params.groups(i);
    
    savename = sprintf('hist_%s.png',group.name);
    relsavename = fullfile(dataloc_params.figdir,savename);
    savename = fullfile(figdir,savename);
     
    if ~DEBUG && plothist < 2 && exist(savename,'file'),
      histfiles{i} = relsavename;
      continue;
    end
    
    if verbose > 1,
      fprintf('Plotting histogram for group %s...\n',group);
    end
    
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
  if verbose > 0,
    fprintf('Outputting summary webpage...\n');
  end
  if exist(statswebpagefile,'file'),
    delete(statswebpagefile);
  end
  MakePerFrameStatsWebpage(statswebpagefile,basename,timeseriesfiles,groups,statfiles,...
    statperflyfiles,stimulus_plotparams.features,stimfiles,stimperflyfiles,...
    videoions,fliesplot,stimtrajfile,videofiles,...
    histfiles,{hist_plot_params.groups.name});
end

%% save info to a mat file

if verbose > 0,
  fprintf('Saving debug info...\n');
end
filename = fullfile(expdir,dataloc_params.plotperframestatsinfomatfilestr);
fprintf(logfid,'Saving debug info to file %s...\n',filename);
ppfsinfo = struct;
ppfsinfo.params = struct;
ppfsinfo.params.stats_plotparams = stats_plotparams;
ppfsinfo.params.stimulus_plotparams = stimulus_plotparams;
ppfsinfo.params.hist_plot_params = hist_plot_params;
ppfsinfo.params.condition_plot_params = condition_plot_params;
if makestimvideos,
    ppfsinfo.params.videoparams = videoparams;
end

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

