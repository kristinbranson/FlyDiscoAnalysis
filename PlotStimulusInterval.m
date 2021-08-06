function [hfig,hfigfly] = PlotStimulusInterval(trx,field,ion,basename,varargin)

defaultcolors = struct('stimbkgd',[0.9087 0.7695 0.7960],...
  'stim',[0.6350 0.0780 0.1840],...
  'off',[0 0 0]);

% parse plotting parameters
[hfig,hfigfly,visible,position,axposition,prestim,poststim,colors,ylim,plotflies,maxnflies] = ...
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
  'maxnflies',inf);

fns = setdiff(fieldnames(defaultcolors),fieldnames(colors));
for i = 1:numel(fns),
  colors.(fns{i}) = defaultcolors.(fns{i});
end

% stimulus periods
ind = trx.getIndicatorLED(1);
non = numel(ind.starton);
assert(ion <= non);

% times to plot
f0 = ind.starton(ion);
f1 = ind.endon(ion);
t0 = trx.movie_timestamps{1}(f0);
t1 = trx.movie_timestamps{1}(f1);
tpre = t0-prestim;
tpost = t1+poststim;
[~,fpre] = min(abs(trx.movie_timestamps{1}-tpre));
[~,fpost] = min(abs(trx.movie_timestamps{1}-tpost));

data = trx.(field);
nflies = numel(data);
ipre = fpre-trx.firstframe+1;
ipost = fpost-trx.firstframe+1;
for fly = 1:nflies,
  data{fly} = padgrab(data{fly},nan,1,1,ipre(fly),ipost(fly));
end
data = cat(1,data{:});
meandata = nanmean(data,1);
ndatafly = sum(~isnan(data),2);
%stddata = nanstd(data,1,1);

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
hax = axes('Position',axposition,'Parent',hfig,'YLim',ylim);

patch(hax,[t0,t0,t1,t1,t0],ylim([1,2,2,1,1]),colors.stimbkgd,'LineStyle','none');
hold(hax,'on');

% hfly = plot(hax,trx.movie_timestamps{1}(fpre:fpost),data,'-','LineWidth',.5);
% flycolors = linspace(.5,.75,nflies)'+[0,0,0];
% for i = 1:nflies,
%   set(hfly(i),'Color',flycolors(i,:));
% end
% plot(hax,[t0,t0],ylim,'-','Color',colors.stim);
% plot(hax,[t1,t1],ylim,'-','Color',colors.stim);
plot(hax,trx.movie_timestamps{1}(fpre:f0-1),meandata(1:f0-fpre),'-','LineWidth',1,'Color',colors.off);
plot(hax,trx.movie_timestamps{1}(f1+1:fpost),meandata(f1-fpre+2:end),'-','LineWidth',1,'Color',colors.off);
plot(hax,trx.movie_timestamps{1}(f0:f1),meandata(f0-fpre+1:f1-fpre+1),'-','LineWidth',1,'Color',colors.stim);
set(hax,'YLim',ylim,'XLim',[tpre,tpost]);

ylabel(hax,field,'Interpreter','none');
xlabel(hax,'Time (s)');
title(hax,{basename,sprintf('Stimulus period %d',ion)},'Interpreter','none');

if plotflies,
  
  nfliesreal = nnz(ndatafly);
  nfliesplot = min(nfliesreal,maxnflies);
  if nfliesplot < nfliesreal,
    [~,order] = sort(ndatafly,1,'descend');
    fliesplot = sort(order(1:nfliesplot));
  else
    fliesplot = find(ndatafly>0);
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
  hax_fly = createsubplots(nfliesplot,1,[[.05,.05];[.05,.01]],hfigfly);

  for flyi = 1:nfliesplot,
    fly = fliesplot(flyi);
  
    patch(hax_fly(flyi),[t0,t0,t1,t1,t0],ylim([1,2,2,1,1]),colors.stimbkgd,'LineStyle','none');
    hold(hax_fly(flyi),'on');
  
    plot(hax_fly(flyi),trx.movie_timestamps{1}(fpre:f0-1),data(fly,1:f0-fpre),'-','LineWidth',1,'Color',colors.off);
    plot(hax_fly(flyi),trx.movie_timestamps{1}(f1+1:fpost),data(fly,f1-fpre+2:end),'-','LineWidth',1,'Color',colors.off);
    plot(hax_fly(flyi),trx.movie_timestamps{1}(f0:f1),data(fly,f0-fpre+1:f1-fpre+1),'-','LineWidth',1,'Color',colors.stim);
    set(hax_fly(flyi),'YLim',ylim,'XLim',[tpre,tpost]);
    
    ylabel(hax_fly(flyi),sprintf('Fly %d',fly),'Interpreter','none');
    if flyi == nfliesplot,
      xlabel(hax_fly(flyi),'Time (s)');
    else
      set(hax_fly(flyi),'XTickLabels',{});
    end
    if flyi == 1,
      title(hax_fly(flyi),{basename,sprintf('%s, period %d',field,ion)},'Interpreter','none');
    end
  end
  
end
