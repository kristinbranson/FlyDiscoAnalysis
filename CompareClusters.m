function [sortedlinesselected,newisin,coeffs,z,hfigs] = CompareClusters(zdatacluster_norm,isin,isout,line_names,shortstatnames)

lambda = .01;

is = isin | isout;
n = nnz(is);
X = zdatacluster_norm(is,:);
X(isnan(X)) = 0;
d = size(X,2);

% LDA to find a direction to order along
y = zeros(1,n);
y(isin(is)) = 1;

mu = nan(2,d);
mu(1,:) = nanmean(X(y==0,:),1);
mu(2,:) = nanmean(X(y==1,:),1);
dmu = diff(mu,1,1);
Sb = dmu'*dmu;
S = nan(d,d,2);
S(:,:,1) = nancov(X(y==0,:),1);
S(:,:,2) = nancov(X(y==1,:),1);
Sw = sum(S,3);
Swreg = (1-lambda)*Sw + eye(d)*lambda;

[coeffs,D] = eigs(Sb,Swreg,1);

proj = X*coeffs;

if nanmean(proj(y==0)) > nanmean(proj(y==1)),
  coeffs = -coeffs;
  proj = -proj;
end

z = abs(dmu);
[~,order] = sort(z,'descend');

sig = nan(2,d);
sig(1,:) = nanstd(X(y==0,:),1,1);
sig(2,:) = nanstd(X(y==1,:),1,1);

%%

x0 = (1:d)-.125;
x1 = (1:d)+.125;

hfigs = [14325,234234];
figure(hfigs(1));
clf;
set(hfigs(1),'Position',[10,10,1600,640]);
hax = axes('Position',[.05,.4,.92,.55]);
plot(x0,X(y==0,order),'.','Color',[.7,.7,.7]);
hold on;
plot(x1,X(y==1,order),'.','Color',[1,.7,.7]);
plot(repmat(x0,2,1),[mu(1,order)+sig(1,order);mu(1,order)-sig(1,order)],'-','Color','k');
plot(repmat(x1,2,1),[mu(2,order)+sig(2,order);mu(2,order)-sig(2,order)],'-','Color',[.6,0,0]);
h = [0,0];
h(1) = plot(x0,mu(1,order),'ko','MarkerFaceColor','k');
h(2) = plot(x1,mu(2,order),'ko','MarkerFaceColor',[.6,0,0]);
legend(h,'out','in');

%lasti = find(znorm(order)<2,1);
lasti = 50;
maxv = max(max(mu(:,order(1:lasti))));
minv = min(min(mu(:,order(1:lasti))));
dv = maxv-minv;
axis([.5,lasti+.5,minv-dv*.05,maxv+dv*.05]);

set(gca,'XTick',1:lasti,'XTickLabel',shortstatnames(order(1:lasti)));
hxtick = rotateticklabel(gca);

% colors = jet(256)*.7;
% maxz = z(order(1));
% minz = z(order(lasti));
% colori = max(1,min(256,floor((log(z)-log(minz))/(log(maxz)-log(minz))*256)+1));
% for i = 1:lasti,
%   set(hxtick(i),'Color',colors(colori(order(i)),:));
% end

box off;

figure(hfigs(2));
clf;
set(hfigs(2),'Position',[10,10,1600,640]);
stati = order(1);
plot(proj(y==0),X(y==0,stati),'k.')
hold on;
plot(proj(y==1),X(y==1,stati),'r.')
legend('in','out');
xlabel('Projection');
h = ylabel(shortstatnames{stati});
set(h,'Interpreter','none');
box off;

input('Zoom in, hit enter, then select the new boundary of the class');
[limproj,~] = ginput(1);
newy = proj > limproj;

isnew = newy == 1 & y' == 0;
if any(isnew),
  fprintf('Lines added:\n')
  fprintf('%s\n',line_names{isnew});
end

isold = newy == 0 & y' == 1;
if any(isold),
  fprintf('Lines removed:\n')
  fprintf('%s\n',line_names{isold});
end

newisin = false(size(isin));
newisin(is) = newy;


[~,lineorder] = sort(proj,'descend');
sortedlines = line_names(lineorder);
sortednewy = newy(lineorder);
sortedlinesselected = sortedlines(sortednewy);
