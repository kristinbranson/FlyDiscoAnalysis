% DebugAssignIdentities

%% set up path
addpath ../misc;
addpath ../filehandling;

%% parameters

loaddata = true;
doflip = false;
trxfilename = 'C:\Data\CourtshipBowl\20120719_shelby_roi07.mat';

% for roii = 1:12,
% 
% trxfilename = ['C:\Data\CourtshipBowl\20120719_shelby_roi',sprintf('%02d.mat',roii)];


T = 10;
nids = 2;
mua = [2;3];
mub = [1;1.5];
sigmaa = ones(nids,1);
sigmab = .5*ones(nids,1);
muarea = mua.*mub*pi;
% too lazy to figure this out :)
sigmaarea = nan(nids,1);
for i = 1:nids,
  tmpa = randn(1,100000)*sigmaa(i) + mua(i);
  tmpb = randn(1,100000)*sigmab(i) + mub(i);
  sigmaarea(i) = std(tmpa.*tmpb*pi);
end
sigmamotion = 1;
vel_dampen = .1;

x0 = zeros(nids,1);
y0 = zeros(nids,1);

%% generate data

if ~loaddata,

a = bsxfun(@plus,bsxfun(@times,randn(nids,T),sigmaa),mua);
b = bsxfun(@plus,bsxfun(@times,randn(nids,T),sigmab),mub);
area = a.*b.*pi;

x = nan(nids,T);
y = nan(nids,T);

t = 1;
x(:,t) = randn(nids,1)*sigmamotion + x0;
y(:,t) = randn(nids,1)*sigmamotion + y0;
t = 2;
x(:,t) = randn(nids,1)*sigmamotion + x(:,t-1);
y(:,t) = randn(nids,1)*sigmamotion + y(:,t-1);

for t = 3:T,
  x(:,t) = randn(nids,1)*sigmamotion + (2-vel_dampen)*x(:,t-1) - (1-vel_dampen)*x(:,t-2);
  y(:,t) = randn(nids,1)*sigmamotion + (2-vel_dampen)*y(:,t-1) - (1-vel_dampen)*y(:,t-2);
end

theta = atan2(diff(y,1,2),diff(x,1,2));
theta = [theta(:,1),theta];

end

%% load data

if loaddata,
  
  load(trxfilename);
  x = cat(1,trx.x);
  y = cat(1,trx.y);
  a = cat(1,trx.a)*2;
  b = cat(1,trx.b)*2;
  theta = cat(1,trx.theta);
  area = a.*b*pi;
  [nids,T] = size(x);
  
end

%% plot

clf;
colors = lines(nids);
h = nan(1,nids);
htrx = nan(1,nids);
for i = 1:nids,
  h(i) = plot(nan,nan,'-','Color',colors(i,:));
  hold on;
  htrx(i) = plot(nan,nan,'.-','Color',colors(i,:));
end
xlim = [min(x(:)),max(x(:))];
ylim = [min(y(:)),max(y(:))];
xlim = xlim + [-1,1]*mean(a(:));
ylim = ylim + [-1,1]*mean(a(:));
axis equal;
axis([xlim,ylim]);


for t = 1:T,
  for i = 1:nids,
    updatefly(h(i),x(i,t),y(i,t),theta(i,t),a(i,t)/2,b(i,t)/2);
    set(htrx(i),'XData',x(i,max(t-50,1):t),'YData',y(i,max(t-50,1):t));
  end
  title(num2str(t));
  pause(.1);
end

%% flip ids

idsperm = nan(nids,T);
xtest = nan(size(x));
ytest = nan(size(y));
atest = nan(size(a));
btest = nan(size(b));
areatest = nan(size(area));
thetatest = nan(size(a));
for t = 1:T,
  if doflip
    idsperm(:,t) = randperm(nids);
  else
    idsperm(:,t) = 1:nids;
  end
  xtest(:,t) = x(idsperm(:,t),t);
  ytest(:,t) = y(idsperm(:,t),t);
  atest(:,t) = a(idsperm(:,t),t);
  btest(:,t) = b(idsperm(:,t),t);
  areatest(:,t) = area(idsperm(:,t),t);
  thetatest(:,t) = theta(idsperm(:,t),t);
end

%% fix identities

[idsfit,muafit,mubfit,muareafit,sigmaafit,sigmabfit,sigmaareafit,sigmamotionfit,...
  cost,niters] = AssignIdentities(xtest,ytest,atest,btest,'vel_dampen',vel_dampen);

all(idsperm(:)==idsfit(:))
% allsame = all(idsfit(1,:)==1);
% allflipped = all(idsfit(1,:)==2);
% fprintf('Roi %d: all same = %d, all flipped = %d\n',roii,allsame,allflipped);
% 
% if ~allsame && ~allflipped,
%   break;
% end


%% plot

xfit = nan(size(x));
yfit = nan(size(y));
afit = nan(size(a));
bfit = nan(size(b));
areafit = nan(size(area));
thetafit = nan(size(theta));
for t = 1:T,
  xfit(:,t) = xtest(idsfit(:,t),t);
  yfit(:,t) = ytest(idsfit(:,t),t);
  afit(:,t) = atest(idsfit(:,t),t);
  bfit(:,t) = btest(idsfit(:,t),t);
  areafit(:,t) = areatest(idsfit(:,t),t);
  thetafit(:,t) = thetatest(idsfit(:,t),t);
end

clf;
%colors = lines(nids);
colors = [0,0,.7;1,0,0];
h = nan(1,nids);
htrx = nan(1,nids);
for i = 1:nids,
  h(i) = plot(nan,nan,'-','Color',colors(i,:));
  hold on;
  htrx(i) = plot(nan,nan,'.-','Color',colors(i,:));
end
xlim = [min(x(:)),max(x(:))];
ylim = [min(y(:)),max(y(:))];
xlim = xlim + [-1,1]*mean(a(:));
ylim = ylim + [-1,1]*mean(a(:));
axis equal;
axis([xlim,ylim]);

for t = 1:T,
  for i = 1:nids,
    updatefly(h(i),xfit(i,t),yfit(i,t),thetafit(i,t),afit(i,t)/2,bfit(i,t)/2);
    set(htrx(i),'XData',xfit(i,max(t-50,1):t),'YData',yfit(i,max(t-50,1):t));
  end
  title(num2str(t));
  pause(.033);
end
