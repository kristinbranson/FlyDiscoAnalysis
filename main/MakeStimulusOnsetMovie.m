function MakeStimulusOnsetMovie(readframe,trx,outfile,ion,fly,ind,varargin)
% ion is the stimulus index caller wants a movie of.

winglw = 1;
bodylw = 2;
trxlw = .5;
colors.off = [0,0,0];
colors.stim = [0.6350    0.0780    0.1840];

[hax,prestim,poststim,minboxwidth,boxborder,...
  interp,fillvalues,DEBUG,downsample] = ...
  myparse(varargin,'hax',[],'prestim',.25,'poststim',1,...
  'minboxwidth',100,'boxborder',10,...
  'interp','bilinear','fillvalues',0,...
  'debug',false,'downsample',1);

if isempty(hax),
  hfig = figure;
  hax = axes('Position',[0,0,1,1],'Parent',hfig);
else
  hfig = get(hax,'Parent');
end

% times to plot
timestamp_from_frame_index = trx.movie_timestamps{1} ;
f0_abs = ind.starton(ion);  % frame index of the first stimulus onset
t0 = timestamp_from_frame_index(f0_abs);
%tpre = t0-prestim;
%tpost = t0+poststim;
%[~,fpre_ideal] = min(abs(timestamp_from_frame_index-tpre));
%[~,fpost_ideal] = min(abs(timestamp_from_frame_index-tpost));
dt = median(diff(timestamp_from_frame_index), 'omitnan') ;
fpre_abs = f0_abs-round(prestim/dt) ;  % frame index in movie.ufmf
fpost_abs = f0_abs+round(poststim/dt) ;  % frame index in movie.ufmf

flyfirstframe = trx(fly).firstframe ;
flyendframe = trx(fly).endframe ;

if flyfirstframe<=fpre_abs && fpre_abs<=flyendframe && flyfirstframe<=fpost_abs && fpost_abs<=flyendframe ,
  % all is well
else
  warning('Not making movie of stimulus index %d onset for fly %d b/c its frame interval is not contained within the frame interval of the trajectory', ...
          ion, fly) ;
  return
end
%fpre = max(trx(fly).firstframe,fpre_ideal);
%fpost = min(trx(fly).endframe,fpost_ideal);

fpre_traj = fpre_abs-flyfirstframe+1 ;  % frame index relative to trajectory
fpost_traj = fpost_abs-flyfirstframe+1 ;  % frame index relative to trajectory  
x = trx(fly).x(fpre_traj:fpost_traj);  % we call the chunk of the trajectory we're cutting out here the "snippet"
y = trx(fly).y(fpre_traj:fpost_traj);
theta = trx(fly).theta(fpre_traj:fpost_traj);
a = trx(fly).a(fpre_traj:fpost_traj);
b = trx(fly).b(fpre_traj:fpost_traj);
xwingl = trx(fly).xwingl(fpre_traj:fpost_traj);
ywingl = trx(fly).ywingl(fpre_traj:fpost_traj);
xwingr = trx(fly).xwingr(fpre_traj:fpost_traj);
ywingr = trx(fly).ywingr(fpre_traj:fpost_traj);
timestamps = timestamp_from_frame_index(fpre_traj:fpost_traj);
f0_traj = f0_abs-flyfirstframe+1;
snippet_frame_count = numel(x) ;

% rotate up
% x0 = x(1);
% y0 = y(1);
theta0 = theta(1);
%nframes = fpost-fpre+1;

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
% x0 = mx-w/2;
% y0 = my-w/2;
% x1 = mx+w/2;
% y1 = my+w/2;

x = x-mx;
y = y-my;
xwingl = xwingl - mx;
ywingl = ywingl - my;
xwingr = xwingr - mx;
ywingr = ywingr - my;

T = [1,0,0;0,1,0;-mx,-my,1];
A = [R,zeros(2,1);zeros(1,2),1]*T;
tform = maketform('affine',double(A));  %#ok<MTFA1> 


%% set up to plot
f0_snippet = f0_traj-fpre_traj+1 ;
cla(hax);
set(hfig,'Units','pixels','Position',[10,10,w,w]);
him = imagesc(zeros(round(w),round(w)),'XData',[-w/2,w/2],'YData',[-w/2,w/2],'Parent',hax,[0,255]);
colormap(hax,'gray');
axis(hax,'image','off');
hold(hax,'on');
htrxpre = plot(hax,x(1:f0_snippet),y(1:f0_snippet),'.-','LineWidth',trxlw,'Color',colors.off*.6+.4);  %#ok<NASGU> 
htrxpost = plot(hax,x(f0_snippet:end),y(f0_snippet:end),'.-','LineWidth',trxlw,'Color',colors.stim*.6+.4);  %#ok<NASGU> 
hbody = drawflyo(0,0,0,1,1,'Parent',hax);
set(hbody,'LineWidth',bodylw,'Color',colors.off);
hwing = plot(hax,[0,0,0],[0,0,0],'-','LineWidth',winglw,'Color',colors.off*.4+.6);
htime = text(-w/2,-w/2,'.0s','Color',colors.off,'Parent',hax,'HorizontalAlignment','left','VerticalAlignment','top');



%% draw frames
if exist(outfile,'file'),
  delete(outfile);
end

snippet_frame_index_from_clip_frame_index = [ fliplr(f0_snippet:-downsample:1) f0_snippet+downsample:downsample:snippet_frame_count ] ;
clip_frame_count = numel(snippet_frame_index_from_clip_frame_index) ;
for clip_frame_index = 1 :  clip_frame_count ,
  sfi = snippet_frame_index_from_clip_frame_index(clip_frame_index) ;  % sfi stands for "snippet frame index"
  movie_frame_index = sfi + fpre_abs - 1 ;  % the first frame of the snippet corresponds to fpre_abs'th frame of the movie
  im = readframe(movie_frame_index) ;
  imsz = size(im);
  imcrop = imtransform(im,tform,interp,'udata',[1,imsz(2)],'vdata',[1,imsz(1)],...
                       'xdata',[-w/2,w/2],'ydata',[-w/2,w/2],'fillvalues',fillvalues);  %#ok<DIMTRNS> 
  set(him,'CData',imcrop);
  updatefly(hbody,x(sfi),y(sfi),theta(sfi),a(sfi),b(sfi));
  if movie_frame_index < f0_abs,
    set(hbody,'Color',colors.off);
    set(htime,'Color',colors.off);
  else
    set(hbody,'Color',colors.stim);
    set(htime,'Color',colors.stim);
  end
  set(hwing,'XData',[xwingl(sfi),x(sfi),xwingr(sfi)],'YData',[ywingl(sfi),y(sfi),ywingr(sfi)]);
  set(htime,'String',sprintf('%+.2fs',timestamps(sfi)-t0));
  
  drawnow;
  
  if ~DEBUG,
    frame = getframe(hax);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if clip_frame_index == 1 ,
      imwrite(imind,cm,outfile,'gif','Loopcount',inf,'DelayTime',0);
    else
      imwrite(imind,cm,outfile,'gif','WriteMode','append','DelayTime',0);
    end
  end
  
end  % for clip_frame_index = 1 :  clip_frame_count
