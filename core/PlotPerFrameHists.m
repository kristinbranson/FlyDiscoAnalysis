function handles = PlotPerFrameHists(field,hist_perframefeatures,...
  histperexp,histperfly,bininfo,hist_plot_params,expname,varargin)

% parse plotting parameters
[hfig,hax,visible,position,axposition,...
  linewidth,stdalpha,linestyle,stdstyle,idx] = ...
  myparse(varargin,'hfig',1,...
  'hax',[],...
  'visible','off',...
  'position',[1 1 1000 500],...
  'axposition',[.1,.1,.85,.8],...
  'linewidth',2,...
  'stdalpha',.2,...
  'linestyle','-',...
  'stdstyle','patch',...
  'idx',[]);

if isfield(hist_plot_params,field),
  plottype = hist_plot_params.(field);
else
  plottype = 'linear';
end

black = ones(1,3)/255;

%% set up figure
if ~isempty(hax) && ishandle(hax),
  hfig = get(hax,'Parent');
else
  if ~ishandle(hfig),
    figure(hfig);
  else
    clf(hfig);
  end
  set(hfig,'Visible',visible,'Position',position,'Renderer','painters');
  hax = axes('Position',axposition,'Parent',hfig,'XColor',black,'YColor',black);
end

% which fields, conditions
if isempty(idx),
  idx = find(strcmp({hist_perframefeatures.field},field));
end
ntypes = numel(idx);
fns = cell(1,ntypes);
flyconditions = cell(1,ntypes);
frameconditions = cell(1,ntypes);
for ii = 1:ntypes,
  i = idx(ii);
  fns{ii} = sprintf('%s_fly%s_frame%s',field,...
    hist_perframefeatures(i).flycondition,...
    hist_perframefeatures(i).framecondition);
  flyconditions{ii} = hist_perframefeatures(i).flycondition;
  frameconditions{ii} = hist_perframefeatures(i).framecondition;
end

% choose colors for each plot
typecolors = nan(ntypes,3);
for i = 1:ntypes,
  fn = sprintf('color_%s_%s',frameconditions{i},flyconditions{i});
  if isfield(hist_plot_params,fn),
    typecolors(i,:) = hist_plot_params.(fn);
  end
end
if any(isnan(typecolors(:,1))),
  newcolors = jet(64)*.7;
  missingcolors = find(isnan(typecolors(:,1)));
  for i = missingcolors',
    [tmp,j] = max(min(dist2(newcolors,typecolors),[],2),[],1);
    typecolors(i,:) = newcolors(j,:);
  end
end


%% select out relevant histogram data
[centers,edges,meanfrac,tmp,stderrfrac,nfliesanalyzed,Zpertype] = ...
  SelectHistData(fns,bininfo,hist_plot_params,plottype,...
  histperfly,histperexp);

%%

hold(hax,'on');

[tmp,order] = sort(Zpertype);
% plot the type stderrs
hstd = nan(1,ntypes);
for typeii = 1:ntypes,
  typei = order(typeii);
  typecolor = typecolors(typei,:);
  if strcmpi(stdstyle,'patch'),
    hstd(typeii) = patch([centers,fliplr(centers)],...
      [meanfrac{typei}+stderrfrac{typei},...
      fliplr(meanfrac{typei}-stderrfrac{typei})],...
      typecolor*stdalpha + (1-stdalpha),'LineStyle','none','parent',hax);
  else
    hstd(typeii) = errorbar(centers,meanfrac{typei},stderrfrac{typei},'.',...
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
  typecolor = typecolors(typei,:);
  htype(typei) = plot(hax,centers,meanfrac{typei},...
    linestyle,'color',typecolor,'linewidth',linewidth);
end

% set axis limits
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

% type strings
if all(strcmp(flyconditions,'any')),
  typestrs = frameconditions;
else
  typestrs = cell(1,ntypes);
  for typei = 1:ntypes,
    typestrs{typei} = sprintf('%s, %s',flyconditions{typei},frameconditions{typei});
  end
end

hleg = legend(htype,typestrs,'Location','Best','Parent',hfig,'Interpreter','none');
hxlabel = xlabel(hax,field,'Interpreter','none');
hylabel = ylabel(hax,'Fraction of frames');
hti = title(hax,{sprintf('Histogram of %s, %s binning',field,plottype),expname},'Interpreter','none');

handles = struct;
handles.htype = htype;
handles.hstd = hstd;
handles.hfig = hfig;
handles.hax = hax;
handles.hleg = hleg;
handles.hxlabel = hxlabel;
handles.hylabel = hylabel;
handles.hti = hti;