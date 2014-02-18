function handles = PlotLinePerFrameStats(statfns,linestats,setstats,allstatsnorm,allnframestotal,metadata,line_name,normcontrolmean,normcontrolstd,basename,varargin)

% parse plotting parameters
[hfig,visible,position,axposition,allfields,...
  allframeconditions,allflyconditions,minnframes,minnexps] = ...
  myparse(varargin,'hfig',1,...
  'visible','off',...
  'position',[1 1 2000 720],...
  'axposition',[.025,.3,.95,.65],...
  'fields',{},...
  'frameconditions',{},...
  'flyconditions',{},...
  'minnframes',200,...
  'minnexps',2);

% set up figure
if ~ishandle(hfig),
  figure(hfig);
end

width = min(position(3),position(3)/100*numel(statfns));
position(3) = width;

clf(hfig);
set(hfig,'Units','pixels','Visible',visible,'Position',position);
hax = axes('Position',axposition,'Parent',hfig);

%% select data

lineidx = find(strcmp(linestats.line_names,line_name));
setidx = find(strcmp({setstats.metadata.line_name},line_name));
nsets = numel(setidx);
expidx = find(strcmp({metadata.line_name},line_name));
nexps = numel(expidx);
[~,exp2setidx] = ismember({metadata(expidx).set},{setstats.metadata(setidx).set});

nfns = numel(statfns);

linemeans = nan(1,nfns);
linestds = nan(1,nfns);
setmeans = nan(nsets,nfns);
expmeans = nan(nexps,nfns);
nframesperexp = zeros(nexps,nfns);
nexpsperset = zeros(nsets,nfns);
nsetsperline = zeros(1,nfns);

for fni = 1:nfns,

  fn = statfns{fni};
  if numel(fn) > 63,
    fn = fn(1:63);
  end
  
  linemeans(fni) = (linestats.normmeans.(fn)(lineidx) - normcontrolmean.(fn)) / normcontrolstd.(fn);
  linestds(fni) = linestats.stds.(fn)(lineidx) / normcontrolstd.(fn);
  setmeans(:,fni) = (setstats.normmeans.(fn)(setidx) - normcontrolmean.(fn)) / normcontrolstd.(fn);
  expmeans(:,fni) = (allstatsnorm.(fn)(expidx) - normcontrolmean.(fn)) / normcontrolstd.(fn);
  nframesperexp(:,fni) = allnframestotal.(fn)(expidx);
  nexpsperset(:,fni) = setstats.nexps.(fn)(setidx);
  nsetsperline(fni) = linestats.nsets.(fn)(lineidx);
  
end

linestderrs = linestds ./ sqrt(nsetsperline);

%% colors

idx = find(~isnan(linemeans));
[~,order] = sort(abs(linemeans(idx)));
[~,order] = sort(order);
colors = zeros(nfns,3)+.5;
colors1 = jet(numel(idx))*.7;
colors(idx,:) = colors1(order,:);


%% plot per-exp and per-set data

plot([-1,nfns+2],[0,0],'k:')

hold(hax,'on');

for i = [3,10],
  plot([-1,nfns+2],[0,0]-i,'k:');
  plot([-1,nfns+2],[0,0]+i,'k:');
end


hexp = nan(nsets,nfns);
hset = nan(nsets,nfns);
maxoff = .25;
offperset = linspace(-maxoff,maxoff,nsets);
expok = nframesperexp >= minnframes;
for fni = 1:nfns,
  
  expy = {};
  sety = [];
  for seti = 1:nsets,
    setycurr = setmeans(seti,fni);
    if isnan(setycurr) || nexpsperset(seti,fni) < minnexps,
      continue;
    end
    expy{end+1} = sort(expmeans(exp2setidx==seti&expok(:,fni)',fni)); %#ok<AGROW>
    sety(end+1) = setycurr; %#ok<AGROW>
  end
  ny = numel(expy);
  if ny == 0,
    continue;
  elseif ny == 1,
    x = fni;
  else
    x = offperset(1:ny);
    x = x - mean(x) + fni;
  end
  for seti = 1:ny,
    hexp(seti,fni) = plot(repmat(x(seti),size(expy{seti})),expy{seti},'.-','Color',(colors(fni,:)+3)/4);
  end
  for seti = 1:ny,
    hset(seti,fni) = plot(x(seti),sety(seti),'*','Color',(colors(fni,:)+1)/2);
  end
  
end

%% plot means and stderrs for the line

hlinestd = nan(1,nfns);
hlinemean = nan(1,nfns);
for fni = 1:nfns,
  hlinestd(fni) = plot(fni+[0,0],linemeans(fni)+linestderrs(fni)*[-1,1],'-','Color',colors(fni,:));
  hlinemean(fni) = plot(fni,linemeans(fni),'ok','MarkerFaceColor',colors(fni,:));
end

%% set axis limits

set(hax,'XLim',[0,nfns+1]);
htmp = findobj(hax,'-property','YData');
ydata = get(htmp,'YData');
maxy = max([ydata{:}]);
miny = min([ydata{:}]);

ylim = [miny-(maxy-miny)*.01,maxy+(maxy-miny)*.01];
set(hax,'YLim',ylim);

%% tick labels

m = regexp(statfns,'^(?<field>.*)_fly(?<flycondition>.*)_frame(?<framecondition>.*)$','names','once');
m = [m{:}];
fields = {m.field};
flyconditions = {m.flycondition};
frameconditions = {m.framecondition};

allflyconditions = unique(flyconditions);
allframeconditions = unique(frameconditions);
if numel(allflyconditions) == 1 && numel(allframeconditions) == 1,
  typestrs = fields;
elseif numel(allflyconditions) == 1,
  typestrs = cell(1,nfns);
  for i = 1:nfns,
    if strcmp(frameconditions{i},'any'),
      typestrs{i} = fields{i};
    else
      typestrs{i} = sprintf('%s,%s',fields{i},frameconditions{i});
    end
  end
elseif numel(allframeconditions) == 1,
  typestrs = cell(1,nfns);
  for i = 1:nfns,
    if strcmp(flyconditions{i},'any'),
      typestrs{i} = fields{i};
    else
      typestrs{i} = sprintf('%s,%s',fields{i},flyconditions{i});
    end
  end
else
  typestrs = cell(1,nfns);
  for i = 1:nfns,
    typestrs{i} = fields{i};
    if ~strcmp(flyconditions{i},'any'),
      typestrs{i} = sprintf('%s,%s',typestrs{i},flyconditions{i});
    end
    if ~strcmp(frameconditions{i},'any'),
      typestrs{i} = sprintf('%s,%s',typestrs{i},frameconditions{i});
    end
  end
end
for i = find(strcmp(fields,'fractime')),
  typestrs{i} = typestrs{i}(numel('fractime')+2:end);
end

set(hax,'XTick',1:nfns,'XTickLabel',typestrs);

%% legend

%xlabel(hax,'Experiment','Interpreter','none');
hy = ylabel(hax,'Control stds','Interpreter','None');
hti = title(hax,sprintf('Mean, stderr of mean for %s',basename),'Interpreter','none');

%% rotate ticks
hx = rotateticklabel(hax,90);
% make sure the ticks don't overlap the x-axis
ex = get(hx(1),'Extent');
y1 = ex(2)+ex(4);
offy = max(0,y1-miny);
for i = 1:numel(hx),
  pos = get(hx(i),'Position');
  pos(2) = pos(2) - offy;
  set(hx(i),'Position',pos,'color',colors(i,:));
end
set(hx,'Interpreter','none');
set(hax,'box','off');
set(hax,'TickLength',[.0025,.0025]);
%% outputs

handles = struct;
handles.hfig = hfig;
handles.hax = hax;
handles.hexp = hexp;
handles.hset = hset;
handles.hlinemean = hlinemean;
handles.hlinestd = hlinestd;
handles.hy = hy;
handles.hti = hti;
handles.hx = hx;
handles.position = position;