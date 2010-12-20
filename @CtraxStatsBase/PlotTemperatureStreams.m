function res = PlotTemperatureStreams(obj,varargin)

[ns,expdirs,condition,hax,hfig,figpos,plotstyleparams,condition2type,legendstyleparams] = ...
  myparse(varargin,'n',[],'expdir',{},'condition',[],...
  'hax',[],'hfig',[],'figpos',[],'plotstyleparams',{},...
  'condition2type',[],'legendstyleparams',{});

% choose experiments to plot
if isempty(ns),
  if isempty(expdirs),
    ns = 1:obj.nexpdirs;
  else
    ns = obj.expdir2n(expdirs);
  end
end

% narrow based on condition
if ~isempty(condition),
  doplot = false(size(ns));
  for i = 1:numel(ns),
    n = ns(i);
    doplot(i) = condition(obj.metadata{n});
  end
end

if ~isempty(condition2type),
  exptype = cell(numel(ns));
  for i = 1:numel(ns),
    n = ns(i);
    exptype{i} = condition2type(obj.metadata{n});
  end
  [exptypes,~,exptypeidx] = unique(exptype);
else
  exptypes = obj.expdir_bases(ns);
  exptypeidx = 1:numel(ns);
end

% get axes to plot in
[hax,hfig] = get_axes(hax,hfig,'figpos',figpos);

ntypes = numel(exptypes);
if ntypes <= 7,
  colors = lines(ntypes);
else
  colors = jet(ntypes)*.5;
end
h = nan(1,numel(ns));
hlegend = nan(1,ntypes);
for i = 1:numel(ns),
  n = ns(i);
  exptypei = exptypeidx(i);
  h(i) = plot(hax,obj.temperaturestreams{n}(:,1)-obj.temperaturestreams{n}(1,1),...
    obj.temperaturestreams{n}(:,2),'.-','color',colors(i,:));
  hlegend(exptypei) = h(i);
  if ~isempty(plotstyleparams),
    set(h(i),plotstyleparams{:});
  end
end

xlabel(hax,'Time (s)');
ylabel(hax,'Temperature (C)');

hlegend = legend(hax,hlegend,exptypes,'Location','BestOutside','interpreter','none');
if ~isempty(legendstyleparams),
  set(hlegend,legendstyleparams{:});
end

res.hax = hax;
res.hfig = hfig;
res.hplot = h;
res.hlegend = hlegend;

