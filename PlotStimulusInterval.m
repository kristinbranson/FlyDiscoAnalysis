function [hfig,hfigfly] = PlotStimulusInterval(trx,field,ions,basename,setname,varargin)

defaultcolors = struct('stimbkgd',[0.927 0.8156 0.8368],...
  'stim',[0.6350 0.0780 0.1840],...
  'off',[0 0 0],...
  'stimstd',[0.8175 0.539 0.592],...
  'offstd',[.7,.7,.7]...
);

% parse plotting parameters
[hfig,hfigfly,visible,position,axposition,prestim,poststim,colors,ylim,plotflies,...
  maxnflies,plotstd,fliesplot,maxdiffintervallength,fontsize] = ...
  myparse(varargin,'hfig',gobjects(1),...
  'hfigfly',gobjects(1),...
  'visible','on',...
  'position',[1 1 800 400],...
  'axposition',[.1,.15,.85,.7],...
  'prestim',1,...
  'poststim',1,...
  'colors',defaultcolors,...
  'ylim',[],...
  'plotflies',true,...
  'maxnflies',inf,...
  'plotstd',true,...
  'fliesplot',[],...
  'maxdiffintervallength',2,...
  'fontsize',8);

fns = setdiff(fieldnames(defaultcolors),fieldnames(colors));
for i = 1:numel(fns),
  colors.(fns{i}) = defaultcolors.(fns{i});
end
nons = numel(ions);

% stimulus periods
ind = trx.getIndicatorLED(1);
nontotal = numel(ind.starton);
assert(all(ions <= nontotal));

% times to plot
f0s = ind.starton(ions);
f1s = ind.endon(ions);
timestamp_from_frame_index = trx.movie_timestamps{1} ;
t0s = timestamp_from_frame_index(f0s);
t1s = timestamp_from_frame_index(f1s);
%tpres = t0s-prestim;
%tposts = t1s+poststim;
%[~,fpres] = min(abs(trx.movie_timestamps{1}'-tpres),[],1);
%[~,fposts] = min(abs(trx.movie_timestamps{1}'-tposts),[],1);
dt = median(diff(timestamp_from_frame_index), 'omitnan') ;
fpres = f0s-round(prestim/dt) ;  % frame index in movie.ufmf
fposts = f1s+round(poststim/dt) ;  % frame index in movie.ufmf
tpres = timestamp_from_frame_index(fpres) ;
tposts = timestamp_from_frame_index(fposts) ;

data = trx.(field);
nflies = numel(data);
ndata = cell(1,nflies);
stddata_perfly = cell(nflies,nons);
for fly = 1:nflies,
  datacurr = cell(nons,1);
  for ioni = 1:nons,
    ipre = fpres(ioni)-trx(fly).firstframe+1;
    ipost = fposts(ioni)-trx(fly).firstframe+1;
    datacurr{ioni} = padgrab(data{fly},nan,1,1,ipre,ipost);
  end
  nscurr = cellfun(@numel,datacurr);
  ncurr = min(nscurr);
  if (max(nscurr) - ncurr > maxdiffintervallength) && (fly == 1),
    warning('Not all interval lengths are the same, using the first %d frames (max length = %d)',ncurr,max(nscurr));
  end
  datacurr = cellfun(@(x) x(1:ncurr), datacurr,'Uni',0);
  datacurr = cat(1,datacurr{:});
  meandatacurr = mean(datacurr,1,'omitnan');
  stddatacurr = std(datacurr,1,1,'omitnan');
  ndatacurr = sum(~isnan(datacurr),1);
  data{fly} = meandatacurr;
  stddata_perfly{fly} = stddatacurr;
  ndata{fly} = ndatacurr;
end
data = cat(1,data{:});
stddata_perfly = cat(1,stddata_perfly{:});
meandata = mean(data,1,'omitnan');
ndata = cat(1,ndata{:});
ndatafly = sum(ndata,2);
if plotstd,
  stddata = std(data,1,1,'omitnan');
end

if isempty(ylim),
  miny = min(data(:));
  % miny = min(miny,min(meandata-stddata));
  maxy = max(data(:));
  % maxy = max(maxy,max(meandata+stddata));  
  ylim = [miny,maxy]+(maxy-miny)*[-.01,.01];
end

% set up figure
if isempty(hfig) || ~ishandle(hfig),
  hfig = figure('Visible',visible);
end

clf(hfig);
set(hfig,'Units','pixels','Visible',visible,'Position',position);
hax = axes('Position',axposition,'Parent',hfig,'YLim',ylim,'FontSize',fontsize);

dfs = f1s-f0s;
[~,baseseti] = min(dfs);
dt = t1s(baseseti)-t0s(baseseti);
f0 = f0s(baseseti);
f1 = f1s(baseseti);
t0 = t0s(baseseti);
fpre = fpres(baseseti);
fpost = fposts(baseseti);
tpre = tpres(baseseti);
tpost = tposts(baseseti);

patch(hax,[0,0,dt,dt,0],ylim([1,2,2,1,1]),colors.stimbkgd,'LineStyle','none');
hold(hax,'on');

if plotstd,
  ts = trx.movie_timestamps{1}(fpre:f0-1) - t0;
  mu = meandata(1:f0-fpre);
  sig = stddata(1:f0-fpre);
  patch(hax,[ts,fliplr(ts)],[mu-sig,fliplr(mu+sig)],colors.offstd,'LineStyle','none');
  ts = trx.movie_timestamps{1}(f1+1:fpost) - t0;
  mu = meandata(f1-fpre+2:end);
  sig = stddata(f1-fpre+2:end);
  patch(hax,[ts,fliplr(ts)],[mu-sig,fliplr(mu+sig)],colors.offstd,'LineStyle','none');
  ts = trx.movie_timestamps{1}(f0:f1) - t0;
  mu = meandata(f0-fpre+1:f1-fpre+1);
  sig = stddata(f0-fpre+1:f1-fpre+1);
  patch(hax,[ts,fliplr(ts)],[mu-sig,fliplr(mu+sig)],colors.stimstd,'LineStyle','none');
end

% hfly = plot(hax,trx.movie_timestamps{1}(fpre:fpost),data,'-','LineWidth',.5);
% flycolors = linspace(.5,.75,nflies)'+[0,0,0];
% for i = 1:nflies,
%   set(hfly(i),'Color',flycolors(i,:));
% end
% plot(hax,[t0,t0],ylim,'-','Color',colors.stim);
% plot(hax,[t1,t1],ylim,'-','Color',colors.stim);
plot(hax,trx.movie_timestamps{1}(fpre:f0-1)-t0,meandata(1:f0-fpre),'-','LineWidth',1,'Color',colors.off);
plot(hax,trx.movie_timestamps{1}(f1+1:fpost)-t0,meandata(f1-fpre+2:end),'-','LineWidth',1,'Color',colors.off);
plot(hax,trx.movie_timestamps{1}(f0:f1)-t0,meandata(f0-fpre+1:f1-fpre+1),'-','LineWidth',1,'Color',colors.stim);
set(hax,'YLim',ylim,'XLim',[tpre-t0,tpost-t0]);

ylabel(hax,field,'Interpreter','none','Fontsize',fontsize);
xlabel(hax,'Time (s)');
title(hax,sprintf('Stimulus period %s',setname),'Interpreter','none');

if plotflies,
  
  if isempty(fliesplot),
    nfliesreal = nnz(ndatafly);
    nfliesplot = min(nfliesreal,maxnflies);
    if nfliesplot < nfliesreal,
      [~,order] = sort(ndatafly,1,'descend');
      fliesplot = sort(order(1:nfliesplot));
    else
      fliesplot = find(ndatafly>0);
    end
  else
    nfliesplot = numel(fliesplot);
  end
  % set up figure
  if isempty(hfigfly) || ~ishandle(hfigfly),
    hfigfly = figure('Visible',visible);
  end
  
  ss = get(0,'ScreenSize');
  position_fly = position;
  position_fly(4) = min(ss(4),position(4)*nfliesplot);
  clf(hfigfly);
  set(hfigfly,'Units','pixels','Visible',visible,'Position',position_fly);
  hax_fly = createsubplots(nfliesplot,1,[[.07,.07];[.025,.01]],hfigfly);

  for flyi = 1:nfliesplot,
    fly = fliesplot(flyi);
  
    patch(hax_fly(flyi),[0,0,dt,dt,0],ylim([1,2,2,1,1]),colors.stimbkgd,'LineStyle','none');
    hold(hax_fly(flyi),'on');
  
    if plotstd,
      ts = trx.movie_timestamps{1}(fpre:f0-1) - t0;
      mu = data(fly,1:f0-fpre);
      sig = stddata_perfly(fly,1:f0-fpre);
      patch(hax_fly(flyi),[ts,fliplr(ts)],[mu-sig,fliplr(mu+sig)],colors.offstd,'LineStyle','none');
      ts = trx.movie_timestamps{1}(f1+1:fpost) - t0;
      mu = data(fly,f1-fpre+2:end);
      sig = stddata_perfly(fly,f1-fpre+2:end);
      patch(hax_fly(flyi),[ts,fliplr(ts)],[mu-sig,fliplr(mu+sig)],colors.offstd,'LineStyle','none');
      ts = trx.movie_timestamps{1}(f0:f1) - t0;
      mu = data(fly,f0-fpre+1:f1-fpre+1);
      sig = stddata_perfly(fly,f0-fpre+1:f1-fpre+1);
      patch(hax_fly(flyi),[ts,fliplr(ts)],[mu-sig,fliplr(mu+sig)],colors.stimstd,'LineStyle','none');
    end

    plot(hax_fly(flyi),trx.movie_timestamps{1}(fpre:f0-1)-t0,data(fly,1:f0-fpre),'-','LineWidth',1,'Color',colors.off);
    plot(hax_fly(flyi),trx.movie_timestamps{1}(f1+1:fpost)-t0,data(fly,f1-fpre+2:end),'-','LineWidth',1,'Color',colors.off);
    plot(hax_fly(flyi),trx.movie_timestamps{1}(f0:f1)-t0,data(fly,f0-fpre+1:f1-fpre+1),'-','LineWidth',1,'Color',colors.stim);
    set(hax_fly(flyi),'YLim',ylim,'XLim',[tpre-t0,tpost-t0]);
    
    ylabel(hax_fly(flyi),sprintf('Fly %d',fly),'Interpreter','none');
    if flyi == nfliesplot,
      xlabel(hax_fly(flyi),'Time (s)');
    else
      set(hax_fly(flyi),'XTickLabels',{});
    end
    if flyi == 1,
      title(hax_fly(flyi),sprintf('%s, period %s',field,setname),'Interpreter','none','FontSize',fontsize);
    end
  end
  
end
