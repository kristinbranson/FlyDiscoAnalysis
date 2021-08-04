function [handles,didplot] = PlotPerFrameHistsGroup(group,hist_perframefeatures,...
  histperexp,histperfly,bins,hist_plot_params,condition_plot_params,expname,varargin)

% parse plotting parameters
[hfig,hax,visible,position,axposition,...
  linewidth,stdalpha,linestyle,stdstyle] = ...
  myparse(varargin,'hfig',gobjects(0),...
  'hax',[],...
  'visible','off',...
  'position',[1 1 1000 500],...
  'axposition',[.1,.1,.85,.8],...
  'linewidth',2,...
  'stdalpha',.2,...
  'linestyle','-',...
  'stdstyle','patch');

plottype = group.plotstyle;

black = ones(1,3)/255;

%% select out relevant histogram data

allfns = fieldnames(histperexp);
allfullfns = cell(size(allfns));
for i = 1:numel(allfns),
  allfullfns{i} = histperexp.(allfns{i}).name;
end

[ism,idx] = ismember(group.features,allfullfns);
assert(all(ism));

% which fields, conditions
ntypes = numel(idx);
fields = cell(1,ntypes);
flyconditions = cell(1,ntypes);
frameconditions = cell(1,ntypes);
fns = allfns(idx);
fullfns = allfullfns(idx);
for i = 1:ntypes,
  m = DecodeStatName(fullfns{i});
  assert(~isempty(m));
  flyconditions{i} = m.flycondition;
  frameconditions{i} = m.framecondition;
  fields{i} = m.field;
end

% choose colors for each plot
[typecolors,typemarkers] = SelectPlotParams(condition_plot_params,frameconditions,flyconditions);

[centers,edges,meanfrac,~,stderrfrac,~,Zpertype] = ...
  SelectHistDataGroup(fns,fields,bins,hist_plot_params,plottype,...
  histperfly,histperexp);

isdata = any([meanfrac{:}]>0);
if ~isdata,
  handles = [];
  didplot = false;
  return;
end


%% set up figure
if ~isempty(hax) && ishandle(hax),
  hfig = get(hax,'Parent');
else
  if isempty(hfig) || ~ishandle(hfig),
    hfig = figure('visible',visible);
  else
    clf(hfig);
  end
  set(hfig,'Visible',visible,'Position',position,'Renderer','painters');
  hax = axes('Position',axposition,'Parent',hfig,'XColor',black,'YColor',black);
end

%%

hold(hax,'on');

[~,order] = sort(Zpertype);
% plot the type stderrs
hstd = nan(1,ntypes);
for typeii = 1:ntypes,
  typei = order(typeii);
  typecolor = typecolors{typei};
  if strcmpi(stdstyle,'patch'),
    hstd(typeii) = patch([centers{typei},fliplr(centers{typei})],...
      [meanfrac{typei}+stderrfrac{typei},...
      fliplr(meanfrac{typei}-stderrfrac{typei})],...
      typecolor*stdalpha + (1-stdalpha),'LineStyle','none','parent',hax);
  else
    hstd(typeii) = errorbar(centers{typei},meanfrac{typei},stderrfrac{typei},'.',...
      'parent',hax,'color',typecolor);
  end
    
end

% put stds on the bottom
if strcmpi(stdstyle,'patch'),
  hchil = get(hax,'children');
  set(hax,'Children',[setdiff(hchil(:),hstd(:));flipud(hstd(:))]);
end

% plot the type means
htype = nan(1,ntypes);
for typeii = 1:ntypes,
  typei = order(typeii);
  typecolor = typecolors{typei};
  htype(typei) = plot(hax,centers{typei},meanfrac{typei},...
    linestyle,'color',typecolor,'linewidth',linewidth);
end

minx = min(cellfun(@(x) x(1),edges));
maxx = max(cellfun(@(x) x(end),edges));

% set axis limits
%if strcmpi(plottype,'log') && centers(1) > 0,
%  set(hax,'XScale','log');
%  set(hax,'XLim',[centers(1),edges(end)]);
%else
set(hax,'XScale','linear');
set(hax,'XLim',[minx,maxx]);
%end

allys = [meanfrac{:}]+[stderrfrac{:}];
maxy = max(allys);
mosty = prctile(allys,99.9);
% make sure there is some non-zero data
if maxy <= 0,
  maxy = 1; mosty = 1;
elseif mosty <= 0,
  mosty = maxy;
end

if ~all(isnan(allys)),
  if mosty / maxy > .9,
    set(hax,'YLim',[-.01*maxy,maxy*1.01]);
  else
    set(hax,'YLim',[-.01*mosty,mosty*1.01]);
  end
end

typestrs = AbbrConditionString(frameconditions,flyconditions);

if all(strcmp(fields,fields{1})),
  
  xtick = get(hax,'XTick');
  xticklabel = cellstr(get(hax,'XTickLabel'))';
  xtick = [xtick,centers{1}];
  xticklabel = [xticklabel,cell(1,numel(centers{1}))];
  [xtick,order] = unique(xtick);
  xticklabel = xticklabel(order);
  set(hax,'XTick',xtick,'XTickLabel',xticklabel);
  id = fields{1};
else
  id = '';  
end
  
% type strings

hleg = legend(htype,typestrs,'Location','Best','Parent',hfig,'Interpreter','none');
hxlabel = xlabel(hax,id,'Interpreter','none');
hylabel = ylabel(hax,'Fraction of frames');
hti = title(hax,{sprintf('%s histogram, %s binning',group.name,plottype),expname},'Interpreter','none');

handles = struct;
handles.htype = htype;
handles.hstd = hstd;
handles.hfig = hfig;
handles.hax = hax;
handles.hleg = hleg;
handles.hxlabel = hxlabel;
handles.hylabel = hylabel;
handles.hti = hti;
handles.position = position;

didplot = true;