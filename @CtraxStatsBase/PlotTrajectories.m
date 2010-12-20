function hax = PlotTrajectories(obj,varargin)

[hax,hfig,figpos,expdirs,flies] = ...
  myparse(varargin,'hax',[],'hfig',[],'figpos',[],...
  'expdirs',obj.expdir_bases,'flies',1:obj.nflies);

[ns,flies] = obj.IntersectFliesExpdirs(flies,expdirs);
nflies = length(flies);

% put all the flies on one figure
n1 = round(sqrt(nflies));
n2 = ceil(nflies/n1);
%n1 = 6;
%n2 = 4;
[hax,hfig] = get_axes(hax,hfig,'figpos',figpos,'axparams',{n1,n2,[[.05,.01];[.05,.05]]});
%hax = reshape(hax,[n1,n2])'; % column major
%hax = hax(:)';
if n1*n2 > nflies,
  delete(hax(nflies+1:end));
  hax = hax(1:nflies);
end

% get bounds for all flies
minx = obj.arena_center_mm(1)-obj.arena_radius_mm;
maxx = obj.arena_center_mm(1)+obj.arena_radius_mm;
miny = obj.arena_center_mm(2)-obj.arena_radius_mm;
maxy = obj.arena_center_mm(2)+obj.arena_radius_mm;
dx = maxx - minx;
dy = maxy - miny;

thetas = linspace(0,2*pi,100);
x_arena = obj.arena_center_mm(1) + obj.arena_radius_mm*cos(thetas);
y_arena = obj.arena_center_mm(2) + obj.arena_radius_mm*sin(thetas);
for i = 1:nflies,
  
  fly = flies(i);
  n = obj.fly2movie(fly);
  
  % plot just the trajectory
  axes(hax(i));
  [bkgdImage,xdata,ydata] = obj.getBkgdImage('n',n);
  image(xdata,ydata,...
    repmat(uint8(bkgdImage),[1,1,3]));
  hold on;
  plot(x_arena,y_arena,'r-');
  hold on;
  plot(obj.trx(fly).x_mm,obj.trx(fly).y_mm,'k.-','markersize',3,'linewidth',.5);
  title(sprintf('Fly %d, %s',fly,obj.expdir_bases{n}),'interpreter','none');
  axis image;
  axis xy;
  axis([minx-.025*dx,maxx+.025*dx,miny-.025*dy,maxy+.025*dy]);
  
  % onlt put an x-axis on the lowest plots
  [c,r] = ind2sub([n1,n2],i);
  if r ~= n1 && r*n2+c <= nflies,
    set(hax(i),'xticklabel',{});
  end
  if c ~= 1,
    set(hax(i),'yticklabel',{});
  end
end

linkaxes(hax);