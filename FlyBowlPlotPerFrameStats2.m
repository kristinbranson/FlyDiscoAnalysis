function FlyBowlPlotPerFrameStats2(expdir,varargin)

special_cases = {'fractime','duration','boutfreq'};

[analysis_protocol,settingsdir,datalocparamsfilestr,visible,controldatadirstr,DEBUG,...
  usedaterange,usepreviousmonth,plotstatsonly,controlstats,controlhist] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'visible','off',...
  'controldatadirstr','current',...
  'debug',false,...
  'usedaterange',true,...
  'usepreviousmonth',false,...
  'plotstatsonly',false,...
  'controlstats',[],...
  'controlhist',[]);

if ischar(plotstatsonly),
  plotstatsonly = str2double(plotstatsonly);
  if isnan(plotstatsonly),
    plotstatsonly = false;
  else
    plotstatsonly = plotstatsonly ~= 0;
  end
end

dateformat = 'yyyymmddTHHMMSS';

%% data locations

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

%% load experiment data

statsmatsavename = fullfile(expdir,dataloc_params.statsperframematfilestr);
load(statsmatsavename,'statsperfly','statsperexp');
if ~plotstatsonly,
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

statsperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statsperframefeaturesfilestr);
stats_perframefeatures = ReadStatsPerFrameFeatures2(statsperframefeaturesfile);

%% load hist params

histperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframefeaturesfilestr);
hist_perframefeatures = ReadHistPerFrameFeatures2(histperframefeaturesfile);
histperframebinsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframebinsfilestr);
load(histperframebinsfile,'bins');
statframeconditionsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.statframeconditionfilestr);
frameconditiondict = ReadParams(statframeconditionsfile);

%% read plotting parameters

histplotparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histplotparamsfilestr);
hist_plot_params = ReadParams(histplotparamsfile);
[tmp,expname] = fileparts(expdir);

%% get control data

if ~(~isempty(controlstats) && ~isempty(controlhist)) && ~isempty(controldatadirstr),
  
  % try to parse date for this experiment
  metadata = parseExpDir(expdir);
  if ~isempty(metadata) && usedaterange,
    dv = datevec(metadata.date,dateformat);
    year = dv(1);
    month = dv(2);
    if usepreviousmonth,
      % get the previous month
      dvend = [year,month,1,0,0,0];
      if month == 1,
        month = 12;
        year = year - 1;
      else
        month = month-1;
      end
      dvstart = [year,month,1,0,0,0];
    else
      % get the current month
      dvstart = [year,month,1,0,0,0];
      if month == 12,
        month = 1;
        year = year + 1;
      else
        month = month+1;
      end
      dvend = [year,month,1,0,0,0];
    end
    daterange = {datestr(dvstart,dateformat),datestr(dvend,dateformat)};
    controldatadir = fullfile(dataloc_params.pBDPGAL4Ustatsdir,sprintf('%sto%s_%s',daterange{:},controldatadirstr));
    if ~exist(controldatadir,'dir'),
      controldatadir = fullfile(dataloc_params.pBDPGAL4Ustatsdir,controldatadirstr);
    end
  else
    controldatadir = fullfile(dataloc_params.pBDPGAL4Ustatsdir,controldatadirstr);
  end
  controlstatsname = fullfile(controldatadir,dataloc_params.statsperframematfilestr);
  controlstats = load(controlstatsname);
  controlhistname = fullfile(controldatadir,dataloc_params.histperframematfilestr);
  if ~plotstatsonly,
    controlhist = load(controlhistname);
  end
  
  % make a soft-link to the control statistics directory
  if isunix && ~DEBUG,
    [~,link] = unix(sprintf('readlink %s',controldatadir));
    if ~isempty(link),
      if link(1) ~= '/',
        realcontroldatadir = fullfile(dataloc_params.pBDPGAL4Ustatsdir,strtrim(link));
      else
        realcontroldatadir = strtrim(link);
      end
    else
      realcontroldatadir = controldatadir;
    end
    if exist(fullfile(figdir,'pBDPGAL4U_stats'),'file'),
      unix(sprintf('rm %s',fullfile(figdir,'pBDPGAL4U_stats')));
    end
    cmd = sprintf('ln -s %s %s',realcontroldatadir,fullfile(figdir,'pBDPGAL4U_stats'));
    unix(cmd);
  end
    
end

%% plot means, stds

[tmp,basename] = fileparts(expdir);
stathandles = PlotPerFrameStats(stats_perframefeatures,statsperfly,statsperexp,controlstats,basename,'visible',visible);
drawnow;  
% plot a second time to get this to work on the cluster
stathandles = PlotPerFrameStats(stats_perframefeatures,statsperfly,statsperexp,controlstats,basename,'visible',visible);
drawnow;  
if ~DEBUG,
  savename = sprintf('stats.png');
  savename = fullfile(figdir,savename);
  if exist(savename,'file'),
    delete(savename);
  end
  set(stathandles.hfig,'Units','pixels','Position',stathandles.position);
  save2png(savename,stathandles.hfig);
end

if plotstatsonly,
  close all;
  return;
end


%% plot histograms

histfields = cell(1,numel(hist_perframefeatures));
histids = cell(1,numel(hist_perframefeatures));
for i = 1:numel(hist_perframefeatures),
  histfields{i} = hist_perframefeatures(i).field;
  if ismember(histfields{i},special_cases),
    histids{i} = [histfields{i},'_',hist_perframefeatures(i).framecondition];
%     frameconditionparams = DecodeConditions(hist_perframefeatures(i).framecondition,frameconditiondict);
%     m = regexp(frameconditionparams(1:2:end),'^[^_]+_(.+)_labels$','once','tokens');
%     tmp = find(~cellfun(@isempty,m),1);
%     histids{i} = [histfields{i},'_',m{tmp}{1}];
  else
    if strcmp(hist_perframefeatures(i).framecondition,'any'),
      histids{i} = histfields{i};
    else
      histids{i} = [histfields{i},'_',hist_perframefeatures(i).framecondition];
    end
  end
end

% histfields = {hist_perframefeatures.field};
% histids = histfields;
% idxdur = find(ismember(histids,special_cases));
% for i = idxdur(:)',
%   frameconditionparams = DecodeConditions(hist_perframefeatures(i).framecondition,frameconditiondict);
%   m = regexp(frameconditionparams(1:2:end),'^[^_]+_(.+)_labels$','once','tokens');
%   tmp = find(~cellfun(@isempty,m),1);
%   histids{i} = [histids{i},'_',m{tmp}{1}];
% end

[histids,tmp,histidx] = unique(histids);
histfields = histfields(tmp);

for i = 1:numel(histids),
  
  field = histfields{i};
  id = histids{i};
  idxcurr = find(histidx == i);
  if strcmp(field,'duration'),
    binfn = id;
  else
    binfn = field;
  end
  
  if ~isempty(controlhist),
    if isfield(controlhist,'meanhistperexp'),
      handles_control = PlotPerFrameHists2Exp(id,field,idxcurr,hist_perframefeatures,...
        controlhist.meanhistperexp,controlhist.histperexp,...
        bins.(binfn),hist_plot_params,expname,...
        'visible',visible,'linestyle',':','stdstyle','errorbar');
    else
      handles_control = PlotPerFrameHists2(id,field,idxcurr,hist_perframefeatures,...
        controlhist.meanhistperfly,controlhist.histperfly,...
        bins.(binfn),hist_plot_params,expname,...
        'visible',visible,'linestyle',':','stdstyle','errorbar');
    end
    hax = handles_control.hax;
  else
    hax = [];
  end
  
  
  
  handles = PlotPerFrameHists2(id,field,idxcurr,hist_perframefeatures,...
    histperexp,histperfly,...
    bins.(binfn),hist_plot_params,expname,...
    'visible',visible,...
    'hax',hax);
  
  if ~isempty(controlhist),
    % fix legend
    s = get(handles.hleg,'String');
    legend([handles_control.htype(1),handles.htype],[{'control'},s],'Parent',handles.hfig,'Interpreter','none');
  end
  
  drawnow;
  if ~DEBUG,
    savename = sprintf('hist_%s.png',id);
    savename = fullfile(figdir,savename);
    if exist(savename,'file'),
      delete(savename);
    end
    set(handles.hfig,'Units','pixels','Position',handles.position);
    save2png(savename,handles.hfig);
  end
  
end

%% plot means, stds
% 
% [tmp,basename] = fileparts(expdir);
% stathandles = PlotPerFrameStats(stats_perframefeatures,statsperfly,statsperexp,controlstats,basename,'visible',visible);
% drawnow;  
% if ~DEBUG,
%   savename = sprintf('stats.png');
%   savename = fullfile(figdir,savename);
%   if exist(savename,'file'),
%     delete(savename);
%   end
%   save2png(savename,stathandles.hfig);
% end

close all;

