function cost = MBLeftRightMatchCost(x,posleft,posright,dataleft,dataright)

alpha = x(1);
beta = x(2);
gamma = x(3);
a = x(4);
b = x(5);
c = x(6);

Rx = [1,0,0
  0,cos(alpha),-sin(alpha)
  0,sin(alpha),cos(alpha)];
Ry = [cos(beta),0,sin(beta)
  0,1,0
  -sin(beta),0,cos(beta)];
Rz = [cos(gamma),-sin(gamma),0
  sin(gamma),cos(gamma),0
  0,0,1];

posleft1 = Rx*Ry*Rz*posleft;
posleft1(1,:)= posleft1(1,:) + a;
posleft1(2,:) = posleft1(2,:) + b;
posleft1(3,:) = posleft1(3,:) + c;

dataright1 = nan(size(posright,2),size(dataleft,2));
for i = 1:size(dataleft,2),  
  interpf = TriScatteredInterp(posleft1(1,:)',posleft1(2,:)',posleft1(3,:)',dataleft(:,i));
  dataright1(:,i) = interpf(posright(1,:)',posright(2,:)',posright(3,:)');
end
dataright1(isnan(dataright1)) = 0;
cost1 = sum(abs(dataright1(:)-dataright(:)));

dataleft1 = nan(size(posleft,2),size(dataleft,2));
for i = 1:size(dataright,2),  
  interpf = TriScatteredInterp(posright(1,:)',posright(2,:)',posright(3,:)',dataright(:,i));
  dataleft1(:,i) = interpf(posleft1(1,:)',posleft1(2,:)',posleft1(3,:)');
end
dataleft1(isnan(dataleft1)) = 0;
cost2 = sum(abs(dataleft1(:)-dataleft(:)));

cost = (cost1+cost2)/2;
fprintf('cost for r = (%f,%f,%f), t = (%f,%f,%f) = %f\n',x,cost);