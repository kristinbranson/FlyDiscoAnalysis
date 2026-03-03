% set up path

[~,computername] = system('hostname');
computername = strtrim(computername);

switch computername,
  
  case 'bransonk-lw1',
    
    addpath E:\Code\JCtrax\misc;
    addpath E:\Code\JCtrax\filehandling;
    addpath('E:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
  
  case 'bransonk-lw2',

    addpath C:\Code\JCtrax\misc;
    addpath C:\Code\JCtrax\filehandling;
    addpath('C:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    
  case 'bransonk-desktop',
    
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
    addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
    addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
    settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';

  otherwise
    
    warning('Unknown computer %s. Paths may not be setup correctly',computername);
    addpath C:\Code\JCtrax\misc;
    addpath C:\Code\JCtrax\filehandling;
    addpath('C:\Code\SAGE\MATLABInterface\Trunk\')
    settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    rootdatadir = 'O:\Olympiad_Screen\fly_bowl\bowl_data';
    
end

%% load data

load('alicehitdata.mat');

%% combine data for the same lines

% parse group names
m = regexp(data_stat.groupsgmr,'^(?<line_name>.*)__Rig(?<rig>.*)__(?<set_datetime>.*)$','names');
% make sure every name matches ok
if any(cellfun(@isempty,m)),
  error('Matching failed for some groups');
end
parsed_group = cell2mat(m);

% find repeats
[data_stat.groupsgmr_unique,~,line_idx] = unique({parsed_group.line_name});

% combine repeats
data_stat.normstatgmr_unique = nan(numel(data_stat.groupsgmr_unique),numel(data_stat.statnames));
for i = 1:numel(data_stat.groupsgmr_unique),
  data_stat.normstatgmr_unique(i,:) = nanmean(data_stat.normstatgmr(line_idx==i,:),1);
end

%% take percentiles to figure out the range of the zscores

r = max(abs(prctile(data_stat.normstatgmr_unique(:),[1,99])));
sigfun = @(x) 2./(1+exp(-x*2./r)) - 1;

%% squash the z-scores

data_stat.squashstatgmr_unique = sigfun(data_stat.normstatgmr_unique);
data_stat.ishit_gmrunique = any(abs(data_stat.normstatgmr_unique) > 1.96,2);
data_stat.maxsquashstatgmr_unique = max(abs(data_stat.normstatgmr_unique),[],2);

%% categories for stats

STATTYPEIDXRANGES = [0,40,43,53,62,88,108,128,154,179];
STATTYPENAMES = {'speed','','turn','straight','dist2close','ddist2close','motion2close','pos2close','arena'};
NSTATTYPES = numel(STATTYPENAMES);
STATTYPECOLORS = lines(NSTATTYPES);

%% do pca

x = data_stat.squashstatgmr_unique;
x(isnan(x)) = 0;
x(isinf(x) & x < 0) = -1;
x(isinf(x) & x > 0) = 1;
[N,D] = size(x);

% pc_coeff(:,d) is the dth principal component
% pc_proj(n,d) is the projection of sample n on component d
% pc_eigs(d) is the eigenvalue for component d
[pc_coeff,pc_proj,pc_eigs,pc_tsquared] = princomp(x);

% compute mean-squared error of pca approx
pc_meansqerr = zeros(1,D);
mux = mean(x,1);
xapprox = repmat(mux,[N,1]);
pc_meansqerr0 = mean(sum((xapprox-x).^2,1));
for d = 1:D,
  xapprox = xapprox + pc_proj(:,d)*pc_coeff(:,d)';
  pc_meansqerr(d) = mean(sum((xapprox-x).^2,1));
end

%% visualize the principal components
hfig = 1;
figure(hfig);
clf;
npcshow = 10;
hax = createsubplots(npcshow+1,1,[.05,.05;.05,.01]);
for i = 1:npcshow,
  hold(hax(i),'on');
  for j = 1:NSTATTYPES,
    x = STATTYPEIDXRANGES(j)+1:STATTYPEIDXRANGES(j+1);
    plot(hax(i),x,pc_coeff(x,i),'.-','color',STATTYPECOLORS(j,:));
  end
  ylabel(hax(i),sprintf('PC %d',i));
end
set([hax(1:end-2),hax(end)],'XTickLabel',{});
linkaxes(hax);
ylim = get(hax(1),'YLim');
y = (ylim(1)*2+ylim(2))/3;
hold(hax(end),'on');
for i = 1:numel(STATTYPENAMES),
  x = [STATTYPEIDXRANGES(i)+1,STATTYPEIDXRANGES(i+1)];
  plot(hax(end),x,y([1,1]),'.-','color',STATTYPECOLORS(i,:));
  text(mean(x),y,STATTYPENAMES{i},'parent',hax(end),...
    'color',STATTYPECOLORS(i,:),'BackgroundColor',get(hfig,'color'),...
    'horizontalalignment','center');
end;
axis(hax(end),'off');
set(hax,'XLim',[0,D+1]);
xlabel(hax(end-1),'Feature');
ylabel(hax(end-1),'Coefficient');
title(hax(1),'Visualization of the first 3 principal components of the z-scores');

%% plot the approximation error
hfig = 2;
figure(hfig);
clf;
npcshow = 10;
plot(0:npcshow,[pc_meansqerr0,pc_meansqerr(1:npcshow)],'k.-');
xlabel('N. principal components');
ylabel('Mean squared error');
title('Error of PCA approximation');

%% visualize the lines projected onto the first  principal components
hfig = 3;
figure(hfig);
clf;
npcshow = 4;
hax = createsubplots(npcshow-1,npcshow-1,[.05,.02;.05,.02]);
hax = reshape(hax,[npcshow-1,npcshow-1]);
for pc2 = 2:npcshow,
  set(hax(pc2-1,pc2-1:npcshow-1),'Visible','off');
  for pc1 = 1:pc2-1,
    scatter(hax(pc2-1,pc1),pc_proj(:,pc1),pc_proj(:,pc2),[],data_stat.maxsquashstatgmr_unique,'.');
    maxx = max(abs(pc_proj(:,pc1)));
    maxy = max(abs(pc_proj(:,pc2)));
    set(hax(pc2-1,pc1),'XLim',1.1*[-maxx,maxx],'YLim',1.1*[-maxy,maxy],'CLim',[0,5]);
  end
end
for pc1 = 1:npcshow-1,
  xlabel(hax(npcshow-1,pc1),sprintf('PC %d',pc1))
end
set(hax(1:end-1,:),'XTickLabel',{});
for pc2 = 2:npcshow,
  ylabel(hax(pc2-1,1),sprintf('PC %d',pc2))
end
set(hax(:,2:end),'YTickLabel',{});

title(hax(1,1),'PCA projection of z-scores');

%% plot all the projections

hfig = 4;
figure(hfig);
clf;
npcshow = 10;

[tmp1,tmp2] = meshgrid(1:npcshow,rand(1,N)*.5-.25);
plot((tmp1+tmp2)',pc_proj(:,1:npcshow)','-');
hold on;
h = plot((tmp1+tmp2)',pc_proj(:,1:npcshow)','.');
set(h(data_stat.ishit_gmrunique,:),'color',[.7,0,0]);
set(h(~data_stat.ishit_gmrunique,:),'color','k');
maxy = max(max(abs(pc_proj(:,1:npcshow))));
set(gca,'XLim',[0,npcshow+1],'YLim',1.1*[-maxy,maxy]);
xlabel('Principal component');
ylabel('Projection');
title('PCA projection of z-scores');

%% try to cluster the hit data

x = data_stat.squashstatgmr_unique(data_stat.ishit_gmrunique,:);
x(isnan(x)) = 0;
x(isinf(x) & x < 0) = -1;
x(isinf(x) & x > 0) = 1;
k = 4;
[kmeans_idx,kmeans_centers] = mykmeans(x,k,'replicates',100);

%% plot the clusters

hfig = 5;
figure(hfig);
clf;
kmeans_colors = lines(k);
hax = createsubplots(k,1,[.05,.05;.05,.01]);

for i = 1:k,
  hold(hax(i),'on');
  idxcurr = kmeans_idx==i;
  for j = 1:NSTATTYPES,
    xj = STATTYPEIDXRANGES(j)+1:STATTYPEIDXRANGES(j+1);
    plot(hax(i),xj,x(idxcurr,xj),'-','color',(1+STATTYPECOLORS(j,:))/2);
    plot(hax(i),xj,kmeans_centers(i,xj),'.-','color',STATTYPECOLORS(j,:),'linewidth',3);
  end
end

set(hax,'XLim',[0,D],'YLim',[-1,1]);

ylim = get(hax(end),'YLim');
y = ylim(1)*.95+ylim(2)*.05;
hold(hax(end),'on');
for i = 1:numel(STATTYPENAMES),
  xj = [STATTYPEIDXRANGES(i)+1,STATTYPEIDXRANGES(i+1)];
  plot(hax(end),xj,y([1,1]),'.-','color',STATTYPECOLORS(i,:));
  text(mean(xj),y,STATTYPENAMES{i},'parent',hax(end),...
    'color',STATTYPECOLORS(i,:),'BackgroundColor',get(hax(end),'color'),...
    'horizontalalignment','center');
end;
set(hax,'XTick',[]);
set(hax(end),'clipping','off');

ylabel(hax(end),'Score');
xlabel(hax(end),'Feature');

%% clustergram

x = data_stat.squashstatgmr_unique(data_stat.ishit_gmrunique,:);
x(isnan(x)) = 0;
x(isinf(x) & x < 0) = -1;
x(isinf(x) & x > 0) = 1;

hfig = 6;
figure(hfig);
clf;
clustergram(x,'RowLabels',{parsed_group(data_stat.ishit_gmrunique).line_name},'ColumnLabels',data_stat.statnames);

%% factor analysis

x = data_stat.squashstatgmr_unique;
x(isnan(x)) = 0;
x(isinf(x) & x < 0) = -1;
x(isinf(x) & x > 0) = 1;
nfactors = 5;

[factor_loadings,factor_vars,factor_rotation,factor_stats,factor_scores] = ...
  factoran(x,nfactors,'rotate','none');

%% visualize the factors
hfig = 7;
figure(hfig);
clf;
hax = createsubplots(nfactors+1,1,[.05,.05;.05,.01]);
for i = 1:nfactors,
  hold(hax(i),'on');
  for j = 1:NSTATTYPES,
    x = STATTYPEIDXRANGES(j)+1:STATTYPEIDXRANGES(j+1);
    plot(hax(i),x,factor_loadings(x,i),'.-','color',STATTYPECOLORS(j,:));
  end
  ylabel(hax(i),sprintf('Factor %d',i));
end
set([hax(1:end-2),hax(end)],'XTickLabel',{});
linkaxes(hax);
ylim = get(hax(1),'YLim');
y = (ylim(1)*2+ylim(2))/3;
hold(hax(end),'on');
for i = 1:numel(STATTYPENAMES),
  x = [STATTYPEIDXRANGES(i)+1,STATTYPEIDXRANGES(i+1)];
  plot(hax(end),x,y([1,1]),'.-','color',STATTYPECOLORS(i,:));
  text(mean(x),y,STATTYPENAMES{i},'parent',hax(end),...
    'color',STATTYPECOLORS(i,:),'BackgroundColor',get(hfig,'color'),...
    'horizontalalignment','center');
end;
axis(hax(end),'off');
set(hax,'XLim',[0,D+1]);
xlabel(hax(end-1),'Feature');
ylabel(hax(end-1),'Loading');
title(hax(1),'Visualization of the factors');