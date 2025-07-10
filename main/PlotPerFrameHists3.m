function handles = PlotPerFrameHists3(id,field,idx,hist_perframefeatures,...
  histcombined,histsingle,bininfo,hist_plot_params,expname,varargin)

% parse plotting parameters
[hfig,hax,visible,position,axposition,...
  linewidth,stdalpha,linestyle,stdstyle,plotstderr,datanames,...
  stattype,meanweighttype,stdweighttype,nframestotal] = ...
  myparse(varargin,'hfig',1,...
  'hax',[],...
  'visible','off',...
  'position',[1 1 1000 500],...
  'axposition',[.1,.1,.85,.8],...
  'linewidth',2,...
  'stdalpha',.2,...
  'linestyle','-',...
  'stdstyle','patch',...
  'plotstderr',true,...
  'datanames',{},...
  'stattype','flymeans',...
  'meanweighttype','nframesfly',...
  'stdweighttype','fracframesfly',...
  'nframestotal',[]);

if isfield(hist_plot_params,id),
  plottype = hist_plot_params.(id);
elseif isfield(hist_plot_params,field),
  plottype = hist_plot_params.(field);
else
  plottype = 'linear';
end

black = ones(1,3)/255;

ndata = numel(histcombined);

%% set up figure
if ~isempty(hax) && ishandle(hax),
  hfig = get(hax,'Parent');
else
  if ~ishandle(hfig),
    figure(hfig);
  else
    clf(hfig);
  end
  set(hfig,'Visible',visible,'Units','pixels','Position',position,'Renderer','painters');
  hax = axes('Position',axposition,'Parent',hfig,'XColor',black,'YColor',black);
end

% which fields, conditions
ntypes = numel(idx);
nplot = ntypes*ndata;

fns = cell(1,ntypes);
flyconditions = cell(1,ntypes);
frameconditions = cell(1,ntypes);
for ii = 1:ntypes,
  i = idx(ii);
  fns{ii} = sprintf('%s_fly%s_frame%s',field,...
    hist_perframefeatures(i).flycondition,...
    hist_perframefeatures(i).framecondition);
  if numel(fns{ii}) > 63,
    fns{ii} = fns{ii}(1:63);
  end
  flyconditions{ii} = hist_perframefeatures(i).flycondition;
  frameconditions{ii} = hist_perframefeatures(i).framecondition;
end

% choose colors for each plot
typecolors = nan(nplot,3);
if ndata == 1,
  for i = 1:ntypes,
    fn = sprintf('color_%s',flyconditions{i});
    if isfield(hist_plot_params,fn),
      typecolors(i,:) = hist_plot_params.(fn);
    end
  end
end
if any(isnan(typecolors(:,1))),
  missingcolors = find(isnan(typecolors(:,1)));
  if numel(missingcolors) == nplot,
    if nplot <= 7,
      typecolors = lines(nplot);
    else
      typecolors = jet(nplot)*.7;
    end
  else
    newcolors = jet(64)*.7;
    for i = missingcolors',
      [~,j] = max(min(dist2(newcolors,typecolors),[],2),[],1);
      typecolors(i,:) = newcolors(j,:);
    end
  end
end

%% select out relevant histogram data
[centers,edges,meanfrac,stdfrac,stderrfrac,~,Zpertype] = ...
  SelectHistData3(fns,bininfo,hist_plot_params,plottype,...
  histsingle,histcombined,'stattype',stattype,...
  'meanweighttype',meanweighttype,'stdweighttype',stdweighttype,...
  'nframestotal',nframestotal);

%%

hold(hax,'on');

% sort so that data types with the highest N are on top
Zpertype1 = sum(Zpertype,1);
[~,typeorder] = sort(Zpertype1);
nplot = ntypes*ndata;
order = reshape(1:nplot,[ndata,ntypes]);
order = order(:,typeorder);
% plot the type stderrs
hstd = nan(ndata,ntypes);
for ii = 1:nplot
  i = order(ii);
  color = typecolors(i,:);
  if plotstderr,
    d = stderrfrac{i};
  else
    d = stdfrac{i};
  end
  if strcmpi(stdstyle,'patch'),
    hstd(ii) = patch([centers,fliplr(centers)],...
      [meanfrac{i}+d,...
      fliplr(meanfrac{i}-d)],...
      color*stdalpha + (1-stdalpha),'LineStyle','none','parent',hax);
  else
    hstd(ii) = errorbar(centers,meanfrac{i},d,'.',...
      'parent',hax,'color',color);
  end
    
end

% % put stds on the bottom
% if strcmpi(stdstyle,'patch'),
%   hchil = get(hax,'children');
%   set(hax,'Children',[setdiff(hchil(:),hstd(:));flipud(hstd(:))]);
% end

% plot the type means
htype = nan(ndata,ntypes);
for ii = 1:nplot,
  i = order(ii);
  color = typecolors(i,:);
  htype(i) = plot(hax,centers,meanfrac{i},...
    linestyle,'color',color,'linewidth',linewidth);
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
if ~all(isnan(allys)),
  if mosty / maxy > .9,
    set(hax,'YLim',[-.01*maxy,maxy*1.01]);
  else
    set(hax,'YLim',[-.01*mosty,mosty*1.01]);
  end
end

xtick = get(hax,'XTick');
xticklabel = cellstr(get(hax,'XTickLabel'))';
xtick = [xtick,centers];
xticklabel = [xticklabel,cell(1,numel(centers))];
[xtick,order] = unique(xtick);
xticklabel = xticklabel(order);
set(hax,'XTick',xtick,'XTickLabel',xticklabel);

% legend
tistr = field;
% type strings
typestrs = cell(1,nplot);
for i = 1:nplot,
  typestrs{i} = '';
end
if ndata > 1,
  for typei = 1:ntypes,
    for datai = 1:ndata,
      i = sub2ind([ndata,ntypes],datai,typei);
      typestrs{i} = [typestrs{i},sprintf('%s, ',datanames{datai})];
    end
  end
elseif ~isempty(datanames),
  tistr = [tistr,', ',datanames{1}];
end
if numel(unique(frameconditions)) > 1,
  for typei = 1:ntypes,
    for datai = 1:ndata,
      i = sub2ind([ndata,ntypes],datai,typei);
      typestrs{i} = [typestrs{i},sprintf('%s, ',frameconditions{typei})];
    end
  end
else
  if ~strcmp(frameconditions{1},'any'),
    tistr = [tistr,', ',frameconditions{1}];
  end
end
if numel(unique(flyconditions)) > 1,
  for typei = 1:ntypes,
    for datai = 1:ndata,
      i = sub2ind([ndata,ntypes],datai,typei);
      typestrs{i} = [typestrs{i},sprintf('%s, ',flyconditions{typei})];
    end
  end
else
  if ~strcmp(flyconditions{1},'any'),
    tistr = [tistr,', ',flyconditions{1}];
  end
end
for i = 1:nplot,
  if ~isempty(typestrs{i}),
    typestrs{i} = typestrs{i}(1:end-2);
  end
end
if nplot > 1,
  hleg = legend(htype,typestrs,'Location','Best','Parent',hfig,'Interpreter','none');
else
  hleg = [];
end
hxlabel = xlabel(hax,field,'Interpreter','none');
hylabel = ylabel(hax,'Fraction of frames');
if plotstderr,
  errstr = 'stderr';
else
  errstr = 'std';
end
if strcmpi(stattype,'flymeans'),
  stattypestr = 'per-fly means';
else
  stattypestr = 'per-exp means';
end
hti = title(hax,{sprintf('Histogram of %s, %s binning, mean and %s of %s',tistr,plottype,errstr,stattypestr),expname},'Interpreter','none');

handles = struct;
handles.htype = htype;
handles.hstd = hstd;
handles.hfig = hfig;
handles.hax = hax;
handles.hleg = hleg;
handles.hxlabel = hxlabel;
handles.hylabel = hylabel;
handles.hti = hti;