function hfig = CompareHistograms(expdirs,fn,varargin)

%% parse arguments
[settingsdir,analysis_protocol,datalocparamsfilestr,...
  hfig,position,axposition,...
  linewidth,stdalpha,colors] = ...
  myparse(varargin,...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'analysis_protocol','current',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'hfig',1,...
  'position',[1 1 1500 500],...
  'axposition',[.1,.1,.85,.8],...
  'linewidth',2,...
  'stdalpha',.2,...
  'colors',zeros(0,3));

nexpdirs = numel(expdirs);

%% file locations
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

%% parse fn into field, fly, frame conditions
m = regexp(fn,'^(?<field>.*)_fly(?<flyconditions>.*)_frame(?<frameconditions>.*)$','names','once');
if isempty(m),
  error('Could not parse %s into field, fly, frame conditions',fn);
end
field = m.field;
%flyconditions = m.flyconditions;
%frameconditions = m.frameconditions;

%% read plotting parameters
histplotparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histplotparamsfilestr);
hist_plot_params = ReadParams(histplotparamsfile);
%histperframefeaturesfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframefeaturesfilestr);
%hist_perframefeatures = ReadHistPerFrameFeatures(histperframefeaturesfile);
histperframebinsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.histperframebinsfilestr);
load(histperframebinsfile,'bins');
bininfo = bins.(field);
if isfield(hist_plot_params,field),
  plottype = hist_plot_params.(field);
else
  plottype = 'linear';
end

%% set up figure
black = ones(1,3)/255;

if ~ishandle(hfig),
  figure(hfig);
else
  clf(hfig);
end
set(hfig,'Position',position);
hax = axes('Position',axposition,'Parent',hfig,'XColor',black,'YColor',black);
hold(hax,'on');

%% choose colors for experiments

colors = cat(1,colors,jet(nexpdirs-size(colors,1))*.7);

%% load data for all experiments

meanfrac = cell(1,nexpdirs);
stdfrac = cell(1,nexpdirs);
stderrfrac = cell(1,nexpdirs);
nfliesanalyzed = nan(1,nexpdirs);

for expdiri = 1:nexpdirs,
  
  expdir = expdirs{expdiri};
  
  % load data
  histmatfile = fullfile(expdir,dataloc_params.histperframematfilestr);
  histdata = load(histmatfile);

  % select out relevant histogram data
  [centers,edges,meanfrac(expdiri),stdfrac(expdiri),...
   stderrfrac(expdiri),nfliesanalyzed(expdiri)] = ...
      SelectHistData({fn},bininfo,hist_plot_params,plottype,...
                     histdata.histperfly,histdata.histperexp);
  
end

%% plot the stderrs
legends = cell(1,nexpdirs);
hstd = nan(1,nexpdirs);
for expdiri = 1:nexpdirs,

  expdir = expdirs{expdiri};
  [~,legends{expdiri}] = fileparts(expdir);

  color = colors(expdiri,:);
  hstd(expdiri) = patch([centers,fliplr(centers)],...
    [meanfrac{expdiri}+stderrfrac{expdiri},...
    fliplr(meanfrac{expdiri}-stderrfrac{expdiri})],...
    color,'facealpha',stdalpha,'LineStyle','none','parent',hax);
end

%% plot the expdir means
hexpdir = nan(1,nexpdirs);
[~,order] = sort(nfliesanalyzed);
for expdirii = 1:nexpdirs,
  expdiri = order(expdirii);
  color = colors(expdiri,:);
  hexpdir(expdiri) = plot(hax,centers,meanfrac{expdiri},...
    '-','color',color,'linewidth',linewidth);
end

%% set axis limits
%if strcmpi(plottype,'log') && centers(1) > 0,
%  set(hax,'XScale','log');
%  set(hax,'XLim',[centers(1),edges(end)]);
%else
set(hax,'XScale','linear');
set(hax,'XLim',[edges(1),edges(end)]);
%end

allys = [meanfrac{:}]+[stderrfrac{:}];
maxy = max(allys);
mosty = prctile(allys,99.9);
if mosty / maxy > .9,
  set(hax,'YLim',[-.01*maxy,maxy*1.01]);
else
  set(hax,'YLim',[-.01*mosty,mosty*1.01]);
end

xtick = get(hax,'XTick');
xticklabel = cellstr(get(hax,'XTickLabel'))';
xtick = [xtick,centers];
xticklabel = [xticklabel,cell(1,numel(centers))];
[xtick,order] = unique(xtick);
xticklabel = xticklabel(order);
set(hax,'XTick',xtick,'XTickLabel',xticklabel);

legend(hexpdir,legends,'Location','EastOutside','Parent',hfig,'Interpreter','none');
xlabel(hax,field,'Interpreter','none');
ylabel(hax,'Fraction of frames');
title(hax,sprintf('Histogram of %s, %s binning',fn,plottype),'Interpreter','none');
