function handles = PlotPerFrameStatsGroup(stats_perframefeatures,statsperfly,statsperexp,basename,varargin)

% parse plotting parameters
[hfig,hfigflies,visible,position,axposition,plotparams,plotflies,fontsize] = ...
  myparse(varargin,'hfig',[],...
  'hfigflies',[],...
  'visible','on',...
  'position',[1 1 800 400],...
  'axposition',[.075,.25,.9,.7],...
  'plotparams',struct,...
  'plotflies',true,...
  'fontsize',8);

% set up figure
if isempty(hfig) || ~ishandle(hfig),
  hfig = figure('Visible',visible);
end

%width = min(position(3),position(3)/20*numel(stats_perframefeatures));
%position(3) = width;

clf(hfig);
set(hfig,'Units','pixels','Visible',visible,'Position',position);
hax = axes('Position',axposition,'Parent',hfig);


if plotflies,
  if isempty(hfigflies) || ~ishandle(hfigflies),
    hfigflies = figure('Visible',visible);
  end
  clf(hfigflies);
  set(hfigflies,'Units','pixels','Visible',visible,'Position',position);
  haxflies = axes('Position',axposition,'Parent',hfigflies);
else
  haxflies = [];
end


%% fns corresponding to fields

shortnames = fieldnames(statsperexp);
fullnames = cell(size(shortnames));
for i = 1:numel(shortnames),
  fullnames{i} = statsperexp.(shortnames{i}).name;
end

fns = stats_perframefeatures;
% i don't think this is working _norm 
m = regexp(fns,'^(.*)_fly(.*)_frame(.*)$','once','tokens');
fields = cellfun(@(x) x{1},m,'UniformOutput',false);
flyconditions = cellfun(@(x) x{2},m,'UniformOutput',false);
frameconditions = cellfun(@(x) x{3},m,'UniformOutput',false);
[ism,idx] = ismember(fns,fullnames);
assert(all(ism));
fns = shortnames(idx);

nfns = numel(fns);
nflies = numel(statsperfly.(fns{1}).mean);
for fni = 2:nfns,
  assert(nflies == numel(statsperfly.(fns{fni}).mean));
end

%% collect data

expmeans = nan(1,nfns);
stdoverflies = nan(1,nfns);
nfliesanalyzed = nan(1,nfns);
flymeans = nan(nfns,nflies);
flyweights = nan(nfns,nflies);

for fni = 1:nfns,

  fn = fns{fni};
  
  expmeans(fni) = statsperexp.(fn).meanmean;
  stdoverflies(fni) = statsperexp.(fn).stdmean;
  
  goodidx = ~isnan(statsperfly.(fn).Z) & ~isnan(statsperfly.(fn).fracframesanalyzed);
  nfliesanalyzed(fni) = sum(statsperfly.(fn).fracframesanalyzed(goodidx));
  assert(~isnan(nfliesanalyzed(fni)));
  
  flymeans(fni,:) = statsperfly.(fn).mean;
  flyweights(fni,:) = statsperfly.(fn).fracframesanalyzed / max(statsperfly.(fn).fracframesanalyzed);
  
end

ndatapts = sum(~isnan(flymeans),1);
isdata = ndatapts > 0;
flyweights(isnan(flymeans)) = 0;

%% colors

[colors,markers] = SelectPlotParams(plotparams,frameconditions,flyconditions);

%% plot experiment mean and standard deviation
hold(hax,'on');
h = gobjects(1,nfns);
hstd = gobjects(1,nfns);
hflies = gobjects(1,nfns);
miny = inf;
maxy = -inf;
nfliesplot = nnz(isdata);

jittermag = .15;
jitter = zeros(1,nflies);
jitter(isdata) = linspace(-jittermag,jittermag,nfliesplot);

for i = 1:nfns,
  fliesplotcurr = ~isnan(flymeans(i,:)) & flyweights(i,:) > 0;
  hflies(i) = scatter(hax,i+jitter(fliesplotcurr),flymeans(i,fliesplotcurr)',8*flyweights(i,fliesplotcurr)',colors{i},markers{i},'filled');
  maxy = max(maxy,max(flymeans(i,isdata)));
  miny = min(miny,min(flymeans(i,isdata)));
  h(i) = plot(hax,i,expmeans(i),markers{i},'color',colors{i},'markerfacecolor',colors{i});
  hstd(i) = plot(hax,[i,i],expmeans(i)+stdoverflies(i)*[-1,1],'-','color',colors{i});
  miny = min(miny,expmeans(i)-stdoverflies(i));
end

%% set axis limits

set(hax,'XLim',[0,nfns+1]);
ylim = [miny-(maxy-miny)*.01,maxy+(maxy-miny)*.01];
set(hax,'YLim',ylim);

%% tick labels

allflyconditions = unique(flyconditions);
allframeconditions = unique(frameconditions);
allfields = unique(fields);
areflyconditions = true;
areframeconditions = true;
arefields = true;
if numel(allflyconditions) == 1,
  areflyconditions = false;
  if ~strcmpi(flyconditions{1},'any'),
    basename = [basename,sprintf(', fly: %s',allflyconditions{1})];
  end
end
if numel(allframeconditions) == 1,
  areframeconditions = false;
  if ~strcmpi(frameconditions{1},'any'),
    basename = [basename,sprintf(', frame: %s',allframeconditions{1})];
  end
end
if numel(allfields) == 1,
  arefields = false;
end
typestrs = repmat({''},[1,nfns]);
for i = 1:nfns,
  isfirst = true;
  if arefields,
    typestrs{i} = fields{i};
    isfirst = false;
  end
  if areflyconditions,
    if ~isfirst,
      typestrs{i} = [typestrs{i},', '];
    end
    isfirst = false;
    typestrs{i} = [typestrs{i},flyconditions{i}];
  end
  if areframeconditions,
    if ~isfirst,
      typestrs{i} = [typestrs{i},', '];
    end
    typestrs{i} = [typestrs{i},frameconditions{i}];
  end
end
  
set(hax,'XTick',1:nfns,'XTickLabel',typestrs,'FontSize',fontsize);
set(hax,'TickLabelInterpreter','none')
%% legend

%xlabel(hax,'Experiment','Interpreter','none');
if ~arefields,
  hy = ylabel(hax,fields{1},'Interpreter','none','FontSize',fontsize);
else 
    hy = [];
end
hti = title(hax,basename,'Interpreter','none','FontSize',fontsize);

set(hax,'box','off');

%% plot per-fly data
if plotflies,

  flycolors = jet(nfliesplot)*.7;
  hold(haxflies,'on');
  hperfly = gobjects(nfliesplot,1);
  sperfly = cell(nfliesplot,1);
  fliesplot = find(isdata);
  for ii = 1:nfliesplot,
    i = fliesplot(ii);
    hperfly(ii) = plot(haxflies,(1:nfns)+jitter(i),flymeans(:,i),'-','Color',flycolors(ii,:));
    idxdatacurr = find(~isnan(flymeans(:,i)) & ~isnan(flyweights(:,i)) & flyweights(:,i)>0);
    scatter(haxflies,idxdatacurr+jitter(i),flymeans(idxdatacurr,i),20*flyweights(idxdatacurr,i),flycolors(ii,:),'o','filled');
    sperfly{ii} = sprintf('Fly %d',i);
  end
  set(haxflies,'XLim',[0,nfns+1]);
  ylim = [miny-(maxy-miny)*.01,maxy+(maxy-miny)*.01];
  set(haxflies,'YLim',ylim);
  set(haxflies,'XTick',1:nfns,'XTickLabel',typestrs);
  set(haxflies,'TickLabelInterpreter','none');
  % hx = rotateticklabel(haxflies,90);
  % % make sure the ticks don't overlap the x-axis
  % ex = get(hx(1),'Extent');
  % y1 = ex(2)+ex(4);
  % offy = max(0,y1-miny);
  % for i = 1:numel(hx),
  %   pos = get(hx(i),'Position');
  %   pos(2) = pos(2) - offy;
  %   set(hx(i),'Position',pos,'color',colors{i});
  % end
  % set(hx,'Interpreter','none');
  set(haxflies,'box','off','FontSize',fontsize);
  if ~arefields,
    hy = ylabel(haxflies,fields{1},'Interpreter','none','FontSize',fontsize);
  else
      hy=[];
  end
  hti = title(haxflies,sprintf('Mean per fly %s',basename),'Interpreter','none','FontSize',fontsize);
  legend(hperfly,sperfly);
  
end

%% outputs

handles = struct;
handles.hfig = hfig;
handles.hax = hax;
handles.h = h;
handles.hy = hy;
handles.hti = hti;
handles.position = position;
handles.hfigflies = hfigflies;
handles.haxflies = haxflies;