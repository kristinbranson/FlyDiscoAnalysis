function MakeStimulusOnsetMovie(readframe,trx,outfile,ion,fly,ind,varargin)

winglw = 1;
bodylw = 2;
trxlw = .5;
defaultcolors.off = [0,0,0];
deafaultcolors.stim = [0.6350    0.0780    0.1840];

[hax,prestim,poststim,minboxwidth,boxborder,...
  interp,fillvalues,DEBUG,downsample,colors] = ...
  myparse(varargin,'hax',[],'prestim',.25,'poststim',1,...
  'minboxwidth',100,'boxborder',10,...
  'interp','bilinear','fillvalues',0,...
  'debug',false,'downsample',1,'colors',defaultcolors);

if isempty(hax),
  hfig = figure;
  hax = axes('Position',[0,0,1,1],'Parent',hfig);
else
  hfig = get(hax,'Parent');
end

% times to plot
f0 = ind.starton(ion);
t0 = trx.movie_timestamps{1}(f0);
tpre = t0-prestim;
tpost = t0+poststim;
[~,fpre] = min(abs(trx.movie_timestamps{1}-tpre));
[~,fpost] = min(abs(trx.movie_timestamps{1}-tpost));
fpre = max(trx(fly).firstframe,fpre);
fpost = min(trx(fly).endframe,fpost);

x = trx(fly).x(fpre:fpost);
y = trx(fly).y(fpre:fpost);
theta = trx(fly).theta(fpre:fpost);
a = trx(fly).a(fpre:fpost);
b = trx(fly).b(fpre:fpost);
xwingl = trx(fly).xwingl(fpre:fpost);
ywingl = trx(fly).ywingl(fpre:fpost);
xwingr = trx(fly).xwingr(fpre:fpost);
ywingr = trx(fly).ywingr(fpre:fpost);
timestamps = trx.movie_timestamps{1}(fpre:fpost);

% rotate up
% x0 = x(1);
% y0 = y(1);
theta0 = theta(f0-fpre+1);
nframes = fpost-fpre+1;

% T = [1,0,0
%   0,1,0
%   -x0,-y0,1];
  
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
x0 = mx-w/2;
y0 = my-w/2;
x1 = mx+w/2;
y1 = my+w/2;

x = x-mx;
y = y-my;
xwingl = xwingl - mx;
ywingl = ywingl - my;
xwingr = xwingr - mx;
ywingr = ywingr - my;

T = [1,0,0;0,1,0;-mx,-my,1];
A = [R,zeros(2,1);zeros(1,2),1]*T;
tform = maketform('affine',double(A));


%% set up to plot

cla(hax);
set(hfig,'Units','pixels','Position',[10,10,w,w]);
him = imagesc(zeros(round(w),round(w)),'XData',[-w/2,w/2],'YData',[-w/2,w/2],'Parent',hax,[0,255]);
colormap(hax,'gray');
axis(hax,'image','off');
hold(hax,'on');
htrxpre = plot(hax,x(1:f0-fpre+1),y(1:f0-fpre+1),'.-','LineWidth',trxlw,'Color',colors.off*.6+.4);
htrxpost = plot(hax,x(f0-fpre+1:end),y(f0-fpre+1:end),'.-','LineWidth',trxlw,'Color',colors.stim*.6+.4);
hbody = drawflyo(0,0,0,1,1,'Parent',hax);
set(hbody,'LineWidth',bodylw,'Color',colors.off);
hwing = plot(hax,[0,0,0],[0,0,0],'-','LineWidth',winglw,'Color',colors.off*.4+.6);
htime = text(-w/2,-w/2,'.0s','Color',colors.off,'Parent',hax,'HorizontalAlignment','left','VerticalAlignment','top');



%% draw frames
if exist(outfile,'file'),
  delete(outfile);
end

frames = [fliplr(f0:-downsample:fpre),f0+downsample:downsample:fpost];
is = frames-fpre+1;
for i = is,
  f = i+fpre-1;
  im = readframe(f);
  imsz = size(im);
  imcrop = imtransform(im,tform,interp,'udata',[1,imsz(2)],'vdata',[1,imsz(1)],...
    'xdata',[-w/2,w/2],'ydata',[-w/2,w/2],'fillvalues',fillvalues);
  set(him,'CData',imcrop);
  updatefly(hbody,x(i),y(i),theta(i),a(i),b(i));
  if f < f0,
    set(hbody,'Color',colors.off);
    set(htime,'Color',colors.off);
  else
    set(hbody,'Color',colors.stim);
    set(htime,'Color',colors.stim);
  end
  set(hwing,'XData',[xwingl(i),x(i),xwingr(i)],'YData',[ywingl(i),y(i),ywingr(i)]);
  set(htime,'String',sprintf('%+.2fs',timestamps(i)-t0));
  
  drawnow;
  
  if ~DEBUG,
    frame = getframe(hax);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == is(1),
      imwrite(imind,cm,outfile,'gif','Loopcount',inf,'DelayTime',0);
    else
      imwrite(imind,cm,outfile,'gif','WriteMode','append','DelayTime',0);
    end
  end
  
end
