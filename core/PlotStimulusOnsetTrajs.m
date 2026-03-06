function [hfig,ions,fliesplot] = PlotStimulusOnsetTrajs(trx,varargin)

trxlw = 1;
winglw = .5;
bodylw = 1;

colors.off = [0,0,0];
colors.stim = [0.6350    0.0780    0.1840];
colors.offtrx = colors.off*.7+.3;
colors.stimtrx = colors.stim*.7+.3;
colors.offbody = colors.off;
colors.stimbody = colors.stim;
colors.offwing = colors.off*.4+.6;
colors.stimwing = colors.stim*.4+.6;

[nfliesplot,fliesplot,ions,nperiodsplot,stimsets,hfig,...
  prestim,poststim,minboxwidth,boxborder,downsample,axwpx,...
  visible,maxnflies] = ...
  myparse_nocheck(varargin,'nfliesplot',[],...
  'fliesplot',[],'ions',[],'nperiodsplot',[],...
  'stimsets',cell(0,2),...
  'hfig',gobjects(0),...
  'prestim',.25,'poststim',1,...
  'minboxwidth',100,'boxborder',10,'downsample',30,'axwpx',200,...
  'visible','on',...
  'maxnflies',inf);

if ~isempty(fliesplot),
  nfliesplot = numel(fliesplot);
end
if isempty(nfliesplot) || nfliesplot < 0,
  nfliesplot = trx.nflies;
end
nfliesplot = min(nfliesplot,maxnflies);

ind = trx.getIndicatorLED(1);
non = numel(ind.starton);

if isempty(ions),
  if isempty(nperiodsplot) || nperiodsplot < 0,
    nperiodsplot = non;
  end
  if isempty(stimsets),
    ions = unique(round(linspace(1,non,nperiodsplot)));
  else
    ions = [];
    for seti = 1:size(stimsets,1),
      idxcurr = unique(round(linspace(1,numel(stimsets{seti,2}),nperiodsplot)));
      ions = [ions,stimsets{seti,2}(idxcurr)]; %#ok<AGROW>
    end
  end
end
nperiodsplot = numel(ions);

if isempty(fliesplot),
  
  if nfliesplot < trx.nflies,
    fliesplot = ChooseFliesPlot(trx,ind,ions,nfliesplot);
    fliesplot = sort(fliesplot);
    nfliesplot = numel(fliesplot);
  else
    fliesplot = 1:trx.nflies;
  end
end

if isempty(hfig) || ~ishandle(hfig),
  hfig = figure('Visible',visible);
else
  clf(hfig);
end

screensz = get(0,'ScreenSize');
figw = axwpx*nperiodsplot/.9;
figh = axwpx*nfliesplot/.9;
r = figw/figh;
maxh = min(screensz(4),screensz(3)/r);
if figh > maxh,
  figh = maxh;
  figw = figh*r;
end
set(hfig,'Units','pixels','Position',[1,1,figw,figh],'Visible',visible);
haxs = createsubplots(nfliesplot,nperiodsplot,[.05,0;.05,0],hfig);
haxs = reshape(haxs,[nfliesplot,nperiodsplot]);

allax = [inf,-inf,inf,-inf];
allts = nan(nperiodsplot,3);
T0 = trx.movie_timestamps{1}(1);
isdata = false(nfliesplot,nperiodsplot);
for ioni = 1:numel(ions),
  ion = ions(ioni);
  % times to plot
  f0 = ind.starton(ion);
  t0 = trx.movie_timestamps{1}(f0);
  tpre = t0-prestim;
  tpost = t0+poststim;
  [~,fpre1] = min(abs(trx.movie_timestamps{1}-tpre));
  [~,fpost1] = min(abs(trx.movie_timestamps{1}-tpost));
  allts(ioni,:) = [tpre,t0,tpost]-T0;

  for fliesploti = 1:numel(fliesplot),
    fly = fliesplot(fliesploti);
    
    if f0 < trx(fly).firstframe || f0 > trx(fly).endframe,
      continue;
    end
    
    hax = haxs(fliesploti,ioni);
    
    fpre = max(trx(fly).firstframe,fpre1);
    fpost = min(trx(fly).endframe,fpost1);
    ipre = fpre-trx(fly).firstframe+1;
    ipost = fpost-trx(fly).firstframe+1;
            
    x = trx(fly).x(ipre:ipost);
    y = trx(fly).y(ipre:ipost);
    theta = trx(fly).theta(ipre:ipost);
    a = trx(fly).a(ipre:ipost);
    b = trx(fly).b(ipre:ipost);
    xwingl = trx(fly).xwingl(ipre:ipost);
    ywingl = trx(fly).ywingl(ipre:ipost);
    xwingr = trx(fly).xwingr(ipre:ipost);
    ywingr = trx(fly).ywingr(ipre:ipost);
    
    % rotate up
    theta0 = theta(f0-fpre+1);
    R = [cos(theta0+pi/2),-sin(theta0+pi/2)
      sin(theta0+pi/2),cos(theta0+pi/2)];
    
    p = [x(:),y(:)]*R;
    x = p(:,1);
    y = p(:,2);
    theta = theta-theta0-pi/2;
    p = [xwingl(:),ywingl(:)]*R;
    xwingl = p(:,1);
    ywingl = p(:,2);
    p = [xwingr(:),ywingr(:)]*R;
    xwingr = p(:,1);
    ywingr = p(:,2);
    maxa = prctile(a,99);
    
    % find a box that fits the trajectories
    minx = min([min(x),min(xwingl),min(xwingr)])-boxborder;
    maxx = max([max(x),max(xwingl),max(xwingr)])+boxborder;
    miny = min([min(y)-2*maxa,min(ywingl),min(ywingr)])-boxborder;
    maxy = max([max(y)+2*maxa,max(ywingl),max(ywingr)])+boxborder;
    mx = (minx+maxx)/2;
    my = (miny+maxy)/2;
    wx = maxx-minx+1;
    wy = maxy-miny+1;
    w = max([wx,wy,minboxwidth]);
    x0 = -w/2;
    y0 = -w/2;
    x1 = +w/2;
    y1 = +w/2;
    
    ax = [x0,x1,y0,y1];
    
    x = x-mx;
    y = y-my;
    xwingl = xwingl - mx;
    ywingl = ywingl - my;
    xwingr = xwingr - mx;
    ywingr = ywingr - my;
    
    % plot
    
    cla(hax);
    hold(hax,'on');
    
    frames = [fliplr(f0-downsample:-downsample:fpre),f0+downsample:downsample:fpost,f0];
    is = frames-fpre+1;
    
    plot(hax,x(1:f0-fpre+1),y(1:f0-fpre+1),'.-','LineWidth',trxlw,'Color',colors.offtrx);
    plot(hax,x(f0-fpre+1:end),y(f0-fpre+1:end),'.-','LineWidth',trxlw,'Color',colors.stimtrx);
    for i = is,
      f = i+fpre-1;
      hwing = plot(hax,[xwingl(i),x(i),xwingr(i)],[ywingl(i),y(i),ywingr(i)],'-','LineWidth',winglw);
      hbody = drawflyo(x(i),y(i),theta(i),a(i),b(i),'Parent',hax);
      set(hbody,'LineWidth',bodylw);
      if f < f0,
        set(hbody,'Color',colors.offbody);
        set(hwing,'Color',colors.offwing);
      else
        set(hbody,'Color',colors.stimbody);
        set(hwing,'Color',colors.stimwing);
      end
    end
    
    axis(hax,'equal');%,'off');
    axis(hax,ax);
    
    allax([1,3]) = min(allax([1,3]),ax([1,3]));
    allax([2,4]) = max(allax([2,4]),ax([2,4]));
    
    isdata(fliesploti,ioni) = true;
    
  end
end
for flyi = 1:nfliesplot,
  fly = fliesplot(flyi);
  ylabel(haxs(flyi,1),sprintf('Fly %d',fly));
end
for ioni = 1:nperiodsplot,
  ion = ions(ioni);
  title(haxs(1,ioni),sprintf('Period %d, %.1fs',ion,allts(ioni,2)),'FontSize',10);
end

set(haxs,'YDir','reverse','XTickLabel',{},'YTickLabel',{});

[flyi,ioni] = ind2sub([nfliesplot,nperiodsplot],find(isdata,1));
plot(haxs(flyi,ioni),.9*allax(2)-[trx.pxpermm,0],.9*allax(3)+[0,0],'k-');
text(haxs(flyi,ioni),.9*allax(2)-trx.pxpermm/2,.9*allax(3)+5,'1 mm','HorizontalAlignment','center','VerticalAlignment','top','FontSize',6);

for i = 1:numel(haxs),
  axis(haxs(i),'equal');
  axis(haxs(i),allax);
end
