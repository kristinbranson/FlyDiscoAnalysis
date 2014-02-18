posleft1 = posleft';
posleft1(1,:) = midline(1)-posleft1(1,:);
posright1 = posright';
posright1(1,:) = posright1(1,:)-midline(1);

mu = mean([posleft1,posright1],2);
posleft1 = bsxfun(@minus,posleft1,mu);
posright1 = bsxfun(@minus,posright1,mu);

clf;
plot3(posleft1(1,:),posleft1(2,:),posleft1(3,:),'k.');
hold on;
plot3(posright1(1,:),posright1(2,:),posright1(3,:),'r.');
axis equal;
grid on;

maxdist = 5;
maxtrans = (maxpos-minpos)/10;
ub = [pi,pi,pi,maxtrans];
lb = -ub;

dxs = -2:2;
dys = -2:2;
dzs = -2:2;
nstates = numel(dxs)*numel(dys)*numel(dzs);
xbest = nan(nstates,6);
parfor i = 1:nstates,
  [xi,yi,zi] = ind2sub([numel(dxs),numel(dys),numel(dzs)],i);
  x0 = [zeros(1,3),dxs(xi),dys(yi),dzs(zi)];
  [xbest(i,:),mincost(i)] = ...
    fmincon(@(x) MBLeftRightMatchCost(x,posleft1,posright1,dataleft,dataright),...
    x0,[],[],[],[],lb,ub,[]);
  fprintf('%s: selected x = %s with cost %f\n',mat2str(x0),mat2str(x),mincost(i));
end

[minmincost,i] = min(mincost);
x = xbest(i,:);
fprintf('rotation = %s, translation = %s, cost = %f\n',mat2str(x(1:3)),mat2str(x(4:6)),minmincost);
% rotation = [-0.00420823056252578 -0.00378717698476799 -0.0117908432055422], translation = [-1.29790142326049 -0.17598211043431 0.352772397189968], cost = 42271.591558

