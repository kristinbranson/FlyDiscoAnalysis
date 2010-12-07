function [rfrac,isonfly] = center_of_rotation2(trk,debug)
if ~exist('debug','var')
  debug = false;
end
N = 100;

% note that if dtheta = 0, this will produce 0

cost = cos(trk.theta);
sint = sin(trk.theta);
dacost = 2*diff(trk.a.*cost);
dbcost = 2*diff(trk.b.*cost);
dasint = 2*diff(trk.a.*sint);
dbsint = 2*diff(trk.b.*sint);
Z = dacost .* dbcost + dbsint.*dasint;
Minv = zeros(2,2,trk.nframes-1);
Minv(1,1,:) = dbcost ./ Z;
Minv(1,2,:) = dbsint ./ Z;
Minv(2,1,:) = -dasint ./ Z;
Minv(2,2,:) = dacost ./ Z;
Minv = reshape(Minv,[4,trk.nframes-1]);
x = [ diff(trk.x) ; diff(trk.y) ];
rfrac = -[Minv(1,:) .* x(1,:) + Minv(3,:) .* x(2,:);...
  Minv(2,:) .* x(1,:) + Minv(4,:) .* x(2,:)];
% when no rotation, set center of rotation to middle of fly by default
rfrac(isnan(rfrac)) = 0;
rfrac0 = rfrac;
isoutofbounds = sum(rfrac.^2,1) > 1;
psi = linspace(0,2*pi,N)';
cospsi = cos(psi);
sinpsi = sin(psi);
idx = find(isoutofbounds);
nout = length(idx);
if nout > 0,
  x1 = repmat(trk.x(idx),[N,1]) + repmat(2*trk.a(idx).*cost(idx),[N,1]).*repmat(cospsi,[1,nout]) - ...
    repmat(2*trk.b(idx).*sint(idx),[N,1]).*repmat(sinpsi,[1,nout]);
  y1 = repmat(trk.y(idx),[N,1]) + repmat(2*trk.a(idx).*sint(idx),[N,1]).*repmat(cospsi,[1,nout]) + ...
    repmat(2*trk.b(idx).*cost(idx),[N,1]).*repmat(sinpsi,[1,nout]);
  x2 = repmat(trk.x(idx+1),[N,1]) + repmat(2*trk.a(idx+1).*cost(idx+1),[N,1]).*repmat(cospsi,[1,nout]) - ...
    repmat(2*trk.b(idx+1).*sint(idx+1),[N,1]).*repmat(sinpsi,[1,nout]);
  y2 = repmat(trk.y(idx+1),[N,1]) + repmat(2*trk.a(idx+1).*sint(idx+1),[N,1]).*repmat(cospsi,[1,nout]) + ...
    repmat(2*trk.b(idx+1).*cost(idx+1),[N,1]).*repmat(sinpsi,[1,nout]);
  d = (x1 - x2).^2 + (y1 - y2).^2;
  [mind,j] = min(d,[],1);
  rfrac(1,idx) = cospsi(j);
  rfrac(2,idx) = sinpsi(j);
end
isonfly = ~isoutofbounds;

if debug,
  ntry = 50;
  psitry = linspace(0,2*pi,ntry+2);
  psitry = psitry(2:end-1);
  [rhotry,psitry] = meshgrid(linspace(0,1,ntry),psitry);
  rtry = [rhotry(:).*cos(psitry(:)),rhotry(:).*sin(psitry(:))]';
  rtry = [rtry,[0;0]];
  colorsplot = 'rg';
  colorsplot2 = 'mc';
  [tmp,order] = sort(-abs(modrange(diff(trk.theta),-pi,pi)));
  for i = order,
    u = [0,0]; v = [0,0];
    uout = [0,0]; vout = [0,0];
    utry = zeros(2,length(rtry)); vtry = zeros(2,length(rtry));
    clf; hold on;
    for j = [i,i+1],
      u(j-i+1) = trk.x(j) + rfrac(1,i)*trk.a(j)*2*cos(trk.theta(j)) - rfrac(2,i)*trk.b(j)*2*sin(trk.theta(j));
      v(j-i+1) = trk.y(j) + rfrac(1,i)*trk.a(j)*2*sin(trk.theta(j)) + rfrac(2,i)*trk.b(j)*2*cos(trk.theta(j));
      uout(j-i+1) = trk.x(j) + rfrac0(1,i)*trk.a(j)*2*cos(trk.theta(j)) - rfrac0(2,i)*trk.b(j)*2*sin(trk.theta(j));
      vout(j-i+1) = trk.y(j) + rfrac0(1,i)*trk.a(j)*2*sin(trk.theta(j)) + rfrac0(2,i)*trk.b(j)*2*cos(trk.theta(j));
      utry(j-i+1,:) = trk.x(j) + rtry(1,:)*trk.a(j)*2*cos(trk.theta(j)) - rtry(2,:)*trk.b(j)*2*sin(trk.theta(j));
      vtry(j-i+1,:) = trk.y(j) + rtry(1,:)*trk.a(j)*2*sin(trk.theta(j)) + rtry(2,:)*trk.b(j)*2*cos(trk.theta(j));
      ellipsedraw(trk.a(j)*2,trk.b(j)*2,trk.x(j),trk.y(j),trk.theta(j),colorsplot(j-i+1));
      plot(trk.x(j)+trk.a(j)*2*cos(trk.theta(j)),trk.y(j)+trk.a(j)*2*sin(trk.theta(j)),'x','color',colorsplot(j-i+1));
      htest(j-i+1) = plot([trk.x(j),u(j-i+1),uout(j-i+1)],[trk.y(j),v(j-i+1),vout(j-i+1)],'o-','color',colorsplot(j-i+1),'markerfacecolor',colorsplot(j-i+1));
    end
    d = (utry(1,:)-utry(2,:)).^2 + (vtry(1,:)-vtry(2,:)).^2;
    [dexp,k] = min(d);
    uexp = utry(:,k);
    vexp = vtry(:,k);
    dexp = sqrt(dexp);
    dtest = sqrt(diff(u).^2 + diff(v).^2);
    dout = sqrt(diff(uout).^2 + diff(vout).^2);
    plot(uexp,vexp,'k.-');
    plot(u,v,'b.-');
    plot(uout,vout,'b.-');
    for j = [i,i+1],
      hexp(j-i+1) = plot([uexp(j-i+1),trk.x(j)],[vexp(j-i+1),trk.y(j)],'o-','color',colorsplot2(j-i+1),'markerfacecolor',colorsplot2(j-i+1));
    end
    ax = nan(1,4);
    ax1 = nan(1,4);
    for j = [i,i+1],
      [ax1(1),ax1(2),ax1(3),ax1(4)] = ellipse_to_bounding_box(trk.x(j),trk.y(j),2*trk.a(j),2*trk.b(j),trk.theta(j));
      ax([1,3]) = min(ax([1,3]),ax1([1,3]));
      ax([2,4]) = max(ax([2,4]),ax1([2,4]));
    end
    ax = ax + [-3,3,-3,3];
    legend([htest,hexp],'analytic1','analytic2','empirical1','empirical2');
    title(sprintf('red = frame %d, green = frame %d, dtest = %f, dout = %f, dexp = %f',i,i+1,dtest,dout,dexp));
    axis equal
    axis(ax)
    input('');
  end
end
