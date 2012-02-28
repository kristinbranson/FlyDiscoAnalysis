function handles = PlotPerFrameHistsComparison(field,flycondition,framecondition,...
  histperexp,histperfly,bininfo,hist_plot_params,datanames,varargin)

% parse plotting parameters
[hfig,hax,visible,position,axposition,...
  linewidth,stdalpha,linestyle,stdstyle] = ...
  myparse(varargin,'hfig',1,...
  'hax',[],...
  'visible','off',...
  'position',[1 1 1000 500],...
  'axposition',[.1,.1,.85,.8],...
  'linewidth',2,...
  'stdalpha',.2,...
  'linestyle','-',...
  'stdstyle','patch');

if isfield(hist_plot_params,field),
  plottype = hist_plot_params.(field);
else
  plottype = 'linear';
end

black = ones(1,3)/255;

ndata = numel(histperexp);

%% set up figure
if ~isempty(hax) && ishandle(hax),
  hfig = get(hax,'Parent');
else
  if ~ishandle(hfig),
    figure(hfig);
  else
    clf(hfig);
  end
  set(hfig,'Visible',visible,'units','pixels','Position',position,'Renderer','painters');
  hax = axes('Position',axposition,'Parent',hfig,'XColor',black,'YColor',black);
end

fn = sprintf('%s_fly%s_frame%s',field,flycondition,framecondition);

% choose colors for each plot
datacolors = lines(ndata);

%% select out relevant histogram data
[centers,edges,meanfrac,~,stderrfrac] = ...
  SelectHistDataComparison(fn,bininfo,hist_plot_params,plottype,...
  histperfly,histperexp);

%%

hold(hax,'on');

% plot the type stderrs
hstd = nan(1,ndata);
for datai = 1:ndata,
  datacolor = datacolors(datai,:);
  if strcmpi(stdstyle,'patch'),
    hstd(datai) = patch([centers,fliplr(centers)],...
      [meanfrac{datai}+stderrfrac{datai},...
      fliplr(meanfrac{datai}-stderrfrac{datai})],...
      datacolor*stdalpha + (1-stdalpha),'LineStyle','none','parent',hax);
  else
    hstd(datai) = errorbar(centers,meanfrac{datai},stderrfrac{datai},'.',...
      'parent',hax,'color',datacolor);
  end
    
end

% put stds on the bottom
if strcmpi(stdstyle,'patch'),
  hchil = get(hax,'children');
  set(hax,'Children',[setdiff(hchil(:),hstd(:));flipud(hstd(:))]);
end

% plot the data means
hdata = nan(1,ndata);
for datai = 1:ndata,
  datacolor = datacolors(datai,:);
  hdata(datai) = plot(hax,centers,meanfrac{datai},...
    linestyle,'color',datacolor,'linewidth',linewidth);
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

hleg = legend(hdata,datanames,'Location','Best','Parent',hfig,'Interpreter','none');
hxlabel = xlabel(hax,field,'Interpreter','none');
hylabel = ylabel(hax,'Fraction of frames');
hti = title(hax,sprintf('Histogram of %s, %s binning',fn,plottype),'Interpreter','none');

handles = struct;
handles.hdata = hdata;
handles.hstd = hstd;
handles.hfig = hfig;
handles.hax = hax;
handles.hleg = hleg;
handles.hxlabel = hxlabel;
handles.hylabel = hylabel;
handles.hti = hti;