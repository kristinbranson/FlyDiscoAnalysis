%% set up path

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
addpath /groups/branson/bransonlab/projects/olympiad/anatomy/fileio;

datafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/CollectedPrimaryPerFrameStatsAndAnatomy20130928.mat';
% imagerydatafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/ImageryData20130725.mat';
% vncdatafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/VNCAnnotations20130905.csv';

anatomydatadir = '/nobackup/branson/AverageAnatomyData20130618';

%% load in data

load(datafile);

%% parameters

maxfraclinesmissingdata = 1;
doanatomyprocessing = false;
statfnset = 'many';
lineset = 'all';
disttransform = 'linearthenlog';
linear2loginflectionpt = 3;

%% choose some statistics

ScriptSetStatsToAnalyze;

%% choose some lines

switch lineset,
  
  case 'all',

    line_names_curr = {...
      '.*'
      };
    
  case 'hasanatomyannotation',

    % choose all lines that have some anatomy data
    
    nmanual = zeros(1,nlines);
    compartments = fieldnames(linestats.int_manual);
    ncompartments = numel(compartments);
    for i = 1:ncompartments,
      fn = compartments{i};
      nmanual = nmanual + ~isnan(linestats.int_manual.(fn));
    end
    fracmanual = nmanual / ncompartments;
    
    lineidxcurr = fracmanual > .5;
    line_names_curr = line_names(lineidxcurr);

end
    

%% filter these out

% stats
statidxcurr = false(1,numel(statfns));
for i = 1:numel(statfnscurr),
  statidxcurr = statidxcurr | ~cellfun(@isempty,regexp(statfns,['^',statfnscurr{i},'$'],'once'));
end
statidxcurr = find(statidxcurr);
statfnscurr = statfns(statidxcurr);
nstatscurr = numel(statidxcurr);

% lines
lineidxcurr = false(numel(line_names),1);
for i = 1:numel(line_names_curr),
  lineidxcurr = lineidxcurr | ~cellfun(@isempty,regexp(line_names,['^',line_names_curr{i},'$'],'once'))';
end
lineidxcurr = find(lineidxcurr);
line_names_curr = line_names(lineidxcurr);
nlinescurr = numel(lineidxcurr);

%% collect this data

datacluster = nan(nlinescurr,nstatscurr);
for ii = 1:nstatscurr,
  i = statidxcurr(ii);
  datacluster(:,ii) = linestats.normmeans.(statfns{i})(lineidxcurr);
end

ncompartments = numel(compartments);
anatdata = nan(nlinescurr,ncompartments);
for i = 1:ncompartments,
  anatdata(:,i) = linestats.int_manual.(compartments{i})(lineidxcurr);
end

%% hand-selected correlation removal

normalizeby = {
%  'fractime_flyany_framestop'
%   'velmag_ctr_flyany_frameany'
%   'dcenter_flyany_frameany'
%   'dist2wall_flyany_frameany'
  };

setiscontrol = strcmp({setstats.metadata.line_name},main_control_line_name);
setdata = nan(nnz(setiscontrol),nstatscurr);
for ii = 1:nstatscurr,
  i = statidxcurr(ii);
  setdata(:,ii) = setstats.normmeans.(statfns{i})(setiscontrol);
end

if isempty(normalizeby),
  datacluster_norm = datacluster;
  setdata_norm = setdata;
else
  [~,idx] = ismember(normalizeby,statfnscurr);
  datacluster_norm = nan(size(datacluster));
  datacluster_norm(:,idx(1)) = datacluster(:,idx(1));
  j = idx(1);
  fprintf('%d: %s, average abs val = %f\n',j,statfnscurr{j},mean(abs(datacluster_norm(:,j))));

  normcoeffs = cell(1,nstatscurr);
  for ii = 2:numel(idx),
    
    is = idx(1:ii-1);
    j = idx(ii);
    [normcoeffs{j},~,datacluster_norm(:,j)] = regress(datacluster(:,j),[datacluster(:,is),ones(nlinescurr,1)]);
    fprintf('%d: %s, average abs val = %f\n',j,statfnscurr{j},mean(abs(datacluster_norm(:,j))));
    
  end

  X = [datacluster(:,idx),ones(nlinescurr,1)];
  for j = setdiff(1:nstatscurr,idx),
    idxgood = ~isnan(datacluster(:,j));
    [normcoeffs{j},~,datacluster_norm(idxgood,j)] = regress(datacluster(idxgood,j),X(idxgood,:));
    fprintf('%d: %s, average abs val = %f\n',j,statfnscurr{j},mean(abs(datacluster_norm(idxgood,j))));
  end
  
  % also normalize the control set stats so that we can z-score

  setdata_norm = nan(size(setdata));

  % first feature
  j = idx(1);
  setdata_norm(:,j) = setdata(:,j);

  % next features
  for ii = 2:numel(idx),
    is = idx(1:ii-1);
    j = idx(ii);
    pred = [setdata(:,is),ones(nnz(setiscontrol),1)]*normcoeffs{j};
    setdata_norm(:,j) = setdata(:,j) - pred;
  end

  % rest of features
  X = [setdata(:,idx),ones(nnz(setiscontrol),1)];
  for j = setdiff(1:nstatscurr,idx),
    pred = X*normcoeffs{j};
    setdata_norm(:,j) = setdata(:,j) - pred;
  end

end
  
% compute mean and standard deviations
munorm = nanmean(setdata_norm,1);
signorm = nanstd(setdata_norm,1,1);

% z-score the line data
zdatacluster_norm = bsxfun(@rdivide,bsxfun(@minus,datacluster_norm,munorm),signorm);
zcontrol_setdata_norm = bsxfun(@rdivide,bsxfun(@minus,setdata_norm,munorm),signorm);

% non-normalized version
mucontrol = nanmean(setdata,1);
sigcontrol = nanstd(setdata,1,1);
zdatacluster = bsxfun(@rdivide,bsxfun(@minus,datacluster,mucontrol),sigcontrol);

% remove nans
zdatacluster_norm_nonan = zdatacluster_norm;
zdatacluster_norm_nonan(isnan(zdatacluster_norm)) = 0;

% 
% % fill in nans with zeros
% zdatacluster_norm_nonan = zdatacluster_norm;
% zdatacluster_norm_nonan(isnan(zdatacluster_norm)) = 0;
% 
% % remove some features so that the X matrix is full-rank
% statidxremove_rank = ismember(statfnscurr,...
%   {'max_wing_angle_flyany_framewingflick'
%   'max_absdwing_angle_flyany_framewingflick'
%   'duration_flyany_framewingflick'
%   'dcenter_flyany_framewingflick'
%   'dell2nose_flyany_framewingflick'
%   'wing_anglel_flyany_frameany'
%   'fractime_flyany_framechase_notwingextension'
%   'wing_angle_imbalance_flyany_frameany'}');

% % here is how I selected features to make X full rank
% tmp = zdatacluster_norm_nonan;
% tmp(:,statidxremove_rank) = [];
% maxrank = rank(tmp);
% tmpnames = shortstatnames(~statidxremove_rank);
% idxcanremove = false(1,size(tmp,2));
% for i = 1:size(tmp,2),
%   [~,~,r] = regress(tmp(:,i),tmp(:,[1:i-1,i+1:size(tmp,2)]));
%   if sum(abs(r)) <= .1,
%     fprintf('%d: %s, regression residual sum = %f\n',i,tmpnames{i},sum(abs(r)));
%   end
% 
%   if rank(tmp(:,[1:i-1,i+1:size(tmp,2)])) == maxrank,
%     fprintf('%d: %s\n',i,tmpnames{i});
%     idxcanremove(i) = true;
%   end
% end

%% remove statistics without enough data

statidxremove = find(sum(isnan(zdatacluster_norm),1) >= nlinescurr*maxfraclinesmissingdata);
statfnscurr0 = statfnscurr;
statidxcurr0 = statidxcurr;
datacluster0 = datacluster;
zdatacluster_norm0 = zdatacluster_norm;
zcontrol_setdata_norm0 = zcontrol_setdata_norm;

statfnscurr(statidxremove) = [];
statidxcurr(statidxremove) = [];
datacluster(:,statidxremove) = [];
zdatacluster_norm(:,statidxremove) = [];
zcontrol_setdata_norm(:,statidxremove) = [];
signorm(statidxremove) = [];
nstatscurr = numel(statidxcurr);

%% create short names for stats and lines for plotting

shortstatnames = statfnscurr;
shortstatnames = regexprep(shortstatnames,'_flyany','');
shortstatnames = regexprep(shortstatnames,'^(.*)_fly(.*)_(.*)','$1_$3_$2');
shortstatnames = regexprep(shortstatnames,'^fractime_frame','fractime_');
shortstatnames = regexprep(shortstatnames,'^duration_frame','duration_');
shortstatnames = regexprep(shortstatnames,'_frameany','');
shortstatnames = regexprep(shortstatnames,'frame','');

shortlinenames = line_names_curr;
shortlinenames = regexprep(shortlinenames,'GMR_','R');
shortlinenames = regexprep(shortlinenames,'_AE_01','');
shortlinenames = regexprep(shortlinenames,'_AD_01','D');

%% non-linear transformation

if strcmpi(disttransform,'linearthenlog'),
  
  % line data
  zdatacluster_transform = nan(size(zdatacluster_norm));
  idx = abs(zdatacluster_norm)<=linear2loginflectionpt;
  zdatacluster_transform(idx) = zdatacluster_norm(idx);
  zdatacluster_transform(~idx) = (linear2loginflectionpt+...
    log(abs(zdatacluster_norm(~idx))-linear2loginflectionpt+1)).*...
    sign(zdatacluster_norm(~idx));

  % control set data
  zcontrol_setdata_transform = nan(size(zcontrol_setdata_norm));
  idx = abs(zcontrol_setdata_norm)<=linear2loginflectionpt;
  zcontrol_setdata_transform(idx) = zcontrol_setdata_norm(idx);
  zcontrol_setdata_transform(~idx) = (linear2loginflectionpt+...
    log(abs(zcontrol_setdata_norm(~idx))-linear2loginflectionpt+1)).*...
    sign(zcontrol_setdata_norm(~idx));
  
elseif strcmpi(disttransform,'log'),

  zdatacluster_transform = log(abs(zdatacluster_norm)).*sign(zdatacluster_norm);

  zcontrol_setdata_transform = log(abs(zcontrol_setdata_norm)).*sign(zcontrol_setdata_norm);

else
  
  zdatacluster_transform = zdatacluster_norm;

  zcontrol_setdata_transform = zcontrol_setdata_norm;

end

%% PCA

X = zdatacluster_transform;
X(isnan(X)) = 0;
mu_pca = mean(X,1);
X = bsxfun(@minus,X,mu_pca);
[coeff,proj_pca,latent,tsquared] = princomp(X);

% find the projection of all the control sets
Xcontrol = zcontrol_setdata_transform;
Xcontrol(isnan(Xcontrol)) = 0;
Xcontrol = bsxfun(@minus,Xcontrol,mu_pca);
projcontrol_pca = Xcontrol*coeff;

% control variance
Scontrol_pca = cov(projcontrol_pca,1);
mucontrol_pca = mean(projcontrol_pca,1);
[a,b,theta] = cov2ell(Scontrol_pca(1:2,1:2));

%% plot the PCA

colorbyfns = {'fractime_flyany_framestop','nflies_close_flyany_frameany','dcenter_flyany_frameany','fractime_flyany_framechase','fractime_flyany_frametouch'};

for fni = 1:numel(colorbyfns),
    
  colorbyfn = colorbyfns{fni};
  
  hfig = fni;
  figure(hfig);
  clf;
  
  x = proj_pca(:,1)-mucontrol_pca(1);
  y = proj_pca(:,2)-mucontrol_pca(2);
  
  maxnstdsplot = 3;
  colors = gray(maxnstdsplot);
  fracdataplot = normcdf([1:maxnstdsplot;-1:-1:-maxnstdsplot],0,1);
  fracdataplot = fracdataplot(1,:)-fracdataplot(2,:);
  fracdataplot(2:end) = fracdataplot(2:end)-fracdataplot(1:end-1);
  
  hcov = nan(1,maxnstdsplot);
  for i = maxnstdsplot:-1:1,
    hcov(i) = ellipsedrawpatch(a*i/2,b*i/2,0,0,theta,1-fracdataplot(i)+zeros(1,3));
    if i == 1,
      hold on;
      hax = gca;
    end
  end
  set(hcov,'LineStyle','--','EdgeColor',[.7,.7,.7]);
  xlim = [min(x),max(x)];
  ylim = [min(y),max(y)];
  dx = diff(xlim);
  dy = diff(ylim);
  xlim = xlim+.05*dx*[-1,1];
  ylim = ylim+.05*dy*[-1,1];
  plot([0,0],ylim,'--','Color',[.7,.7,.7]);
  plot(xlim,[0,0],'--','Color',[.7,.7,.7]);
  
  %plot(projcontrol_pca(:,1),projcontrol_pca(:,2),'kx');
  %hold on;
  
  stati = find(strcmp(statfnscurr,colorbyfn));
  ncolors = 256;
  colors = jet(ncolors)*.9;
  minz = min(zdatacluster_transform(:,stati));
  maxz = max(zdatacluster_transform(:,stati));
  colori = max(1,ceil((zdatacluster_transform(:,stati)-minz)/(maxz-minz)*ncolors));
  hpts = nan(1,nlinescurr);
  for i = 1:nlinescurr,
    hpts(i) = plot(x(i),y(i),'.','Color',colors(colori(i),:));
  end
  %hpts = scatter(x,y,[],zdatacluster_transform(:,stati),'.');
  axis equal;
  axis([xlim,ylim]);
  
  xlabel('PC 1');
  ylabel('PC 2');
  
  set(hax,'UserData',struct('selectedpts',[]));
  title(colorbyfn,'Interpreter','none');
  
  hdata_b1b2a1 = SetUpButtonDown_ReturnPointIndex(hax,x,y,{@ButtonDownFcn_SelectAndShowLineInfo,line_names_curr,hpts});

end
  
%DrawBehaviorVisualization(stats,colors,statlims,x1,x2,y1,y2,'hax',hax)


%% sparse pca

X = zdatacluster_transform;
X(isnan(X)) = 0;
mu_spca = mean(X,1);
X = bsxfun(@minus,X,mu_spca);
[coeffs_spca SD L D paths] = spca(X,[],2,0,-10,100,[],true);

for i = 1:size(coeffs_spca,2);
  idx = find(abs(coeffs_spca(:,i)) > 0);
  fprintf('\nPC %d:\n',i);
  [~,order] = sort(abs(coeffs_spca(idx,i)),1,'descend');
  for j = idx(order)',
    fprintf('%s: %f\n',shortstatnames{j},coeffs_spca(j,i));
  end
end

proj_spca = X*coeffs_spca;

% find the projection of all the control sets
Xcontrol = zcontrol_setdata_transform;
Xcontrol(isnan(Xcontrol)) = 0;
Xcontrol = bsxfun(@minus,Xcontrol,mu_spca);
projcontrol_spca = Xcontrol*coeffs_spca;

% control variance
Scontrol_spca = cov(projcontrol_spca,1);
mucontrol_spca = mean(projcontrol_spca,1);
[a,b,theta] = cov2ell(Scontrol_spca(1:2,1:2));


%% plot the SPCA

colorbyfns = {'absdv_cor_flyany_framenearfly','dnose2ell_flyany_frameany'};

for fni = 1:numel(colorbyfns),
    
  colorbyfn = colorbyfns{fni};
  
  hfig = fni + 10;
  figure(hfig);
  clf;
  
  x = proj_spca(:,1)-mucontrol_spca(1);
  y = proj_spca(:,2)-mucontrol_spca(2);
  
  maxnstdsplot = 3;
  colors = gray(maxnstdsplot);
  fracdataplot = normcdf([1:maxnstdsplot;-1:-1:-maxnstdsplot],0,1);
  fracdataplot = fracdataplot(1,:)-fracdataplot(2,:);
  fracdataplot(2:end) = fracdataplot(2:end)-fracdataplot(1:end-1);
  
  hcov = nan(1,maxnstdsplot);
  for i = maxnstdsplot:-1:1,
    hcov(i) = ellipsedrawpatch(a*i/2,b*i/2,0,0,theta,1-fracdataplot(i)+zeros(1,3));
    if i == 1,
      hold on;
      hax = gca;
    end
  end
  set(hcov,'LineStyle','--','EdgeColor',[.7,.7,.7]);
  xlim = [min(x),max(x)];
  ylim = [min(y),max(y)];
  dx = diff(xlim);
  dy = diff(ylim);
  xlim = xlim+.05*dx*[-1,1];
  ylim = ylim+.05*dy*[-1,1];
  plot([0,0],ylim,'--','Color',[.7,.7,.7]);
  plot(xlim,[0,0],'--','Color',[.7,.7,.7]);
  
  %plot(projcontrol_pca(:,1),projcontrol_pca(:,2),'kx');
  %hold on;
  
  stati = find(strcmp(statfnscurr,colorbyfn));
  ncolors = 256;
  colors = jet(ncolors)*.9;
  minz = min(zdatacluster_transform(:,stati));
  maxz = max(zdatacluster_transform(:,stati));
  colori = max(1,ceil((zdatacluster_transform(:,stati)-minz)/(maxz-minz)*ncolors));
  hpts = nan(1,nlinescurr);
  for i = 1:nlinescurr,
    hpts(i) = plot(x(i),y(i),'.','Color',colors(colori(i),:));
  end
  %hpts = scatter(x,y,[],zdatacluster_transform(:,stati),'.');
  axis equal;
  axis([xlim,ylim]);
  
  xlabel('PC 1');
  ylabel('PC 2');
  
  set(hax,'UserData',struct('selectedpts',[]));
  title(colorbyfn,'Interpreter','none');
  
  hdata_b1b2a1 = SetUpButtonDown_ReturnPointIndex(hax,x,y,{@ButtonDownFcn_SelectAndShowLineInfo,line_names_curr,hpts});

end
  
%DrawBehaviorVisualization(stats,colors,statlims,x1,x2,y1,y2,'hax',hax)

%% t-sne projection


X = [zdatacluster_transform;zcontrol_setdata_transform];
X(isnan(X)) = 0;
mu_tsne = mean(X,1);
X = bsxfun(@minus,X,mu_tsne);

stati = find(strcmp(statfnscurr,'frame_flyany_framestop'));
ydata = tsne(X,zdatacluster_transform(:,stati));
proj_tsne = ydata(1:nlines,:);
projcontrol_tsne = ydata(nlines+1:end,:);

% approximate projection for control sets using nearest neighbor
% Xcontrol = zcontrol_setdata_transform;
% Xcontrol(isnan(Xcontrol)) = 0;
% Xcontrol = bsxfun(@minus,Xcontrol,mu_tsne);
% D = dist2(X,Xcontrol);
% [mind,mini] = min(D,[],1);
% 
% projcontrol_tsne = ydata(mini,:);

% control variance
Scontrol_tsne = cov(projcontrol_tsne,1);
mucontrol_tsne = mean(projcontrol_tsne,1);
[a,b,theta] = cov2ell(Scontrol_tsne(1:2,1:2));


%% plot the t-sne projection

colorbyfns = {'fractime_flyany_framestop','nflies_close_flyany_frameany','dcenter_flyany_frameany','fractime_flyany_framechase','fractime_flyany_frametouch','fractime_flyany_framejump','fractime_flyany_framewinggrooming','fractime_flyany_framepivotcenter'};

for fni = 1:numel(colorbyfns),
    
  colorbyfn = colorbyfns{fni};
  
  hfig = fni + 10;
  figure(hfig);
  clf;
  set(hfig,'Units','pixels','Position',[376 1028 382 258]);
  
  x = ydata(:,1) - mucontrol_tsne(1);
  y = ydata(:,2) - mucontrol_tsne(2);
  xcontrol = projcontrol_tsne(:,1)-mucontrol_tsne(1);
  ycontrol = projcontrol_tsne(:,2)-mucontrol_tsne(2);
  
  plot(xcontrol,ycontrol,'k.');
  hold on;
%   
%   maxnstdsplot = 3;
%   colors = gray(maxnstdsplot);
%   fracdataplot = normcdf([1:maxnstdsplot;-1:-1:-maxnstdsplot],0,1);
%   fracdataplot = fracdataplot(1,:)-fracdataplot(2,:);
%   fracdataplot(2:end) = fracdataplot(2:end)-fracdataplot(1:end-1);
%   
%   hcov = nan(1,maxnstdsplot);
%   for i = maxnstdsplot:-1:1,
%     hcov(i) = ellipsedrawpatch(a*i/2,b*i/2,0,0,theta,1-fracdataplot(i)+zeros(1,3));
%     if i == 1,
%       hold on;
%       hax = gca;
%     end
%   end
%   set(hcov,'LineStyle','--','EdgeColor',[.7,.7,.7]);
  xlim = [min(x),max(x)];
  ylim = [min(y),max(y)];
  dx = diff(xlim);
  dy = diff(ylim);
  xlim = xlim+.05*dx*[-1,1];
  ylim = ylim+.05*dy*[-1,1];
  plot([0,0],ylim,'--','Color',[.7,.7,.7]);
  hold on;
  hax = gca;
  plot(xlim,[0,0],'--','Color',[.7,.7,.7]);
  
  %plot(projcontrol_pca(:,1),projcontrol_pca(:,2),'kx');
  %hold on;

  stati = find(strcmp(statfnscurr,colorbyfn));
  
  if false
    ncolors = 256;
    colors = jet(ncolors)*.9;
    minz = min(zdatacluster_transform(:,stati));
    maxz = max(zdatacluster_transform(:,stati));
    colori = max(1,ceil((zdatacluster_transform(:,stati)-minz)/(maxz-minz)*ncolors));
    hpts = nan(1,nlinescurr);
    for i = 1:nlinescurr,
      hpts(i) = plot(x(i),y(i),'.','Color',colors(colori(i),:));
    end
  else
    idx = ismember(line_names,linenames_jump_avoid);
    s = 24+zeros(1,nlinescurr);
    s(idx) = 48;
    hpts = scatter(x(1:nlinescurr),y(1:nlinescurr),[],zdatacluster_norm(:,stati),'.');
    %set(hpts,'MarkerEdgeColor','k')
    box off;
    %set(gca,'Color',[.5,.5,.5]);
    colormap(jet(256)*.75);
    set(gca,'CLim',[-10,10]);
    colorbar;
%     for i = find(idx),
%       text(x(i),y(i),shortlinenames{i},'HorizontalAlignment','left','VerticalAlignment','middle');
%     end
  end
  axis square;
  axis([xlim,ylim]);
  
%   xlabel('PC 1');
%   ylabel('PC 2');
  
  set(hax,'UserData',struct('selectedpts',[]));
  title(colorbyfn,'Interpreter','none');
  
  drawnow;
  
  SaveFigLotsOfWays(hfig,sprintf('projectionfigs/TSNEColoredBy_%s_20130928',colorbyfn));
  %hdata_b1b2a1 = SetUpButtonDown_ReturnPointIndex(hax,x,y,{@ButtonDownFcn_SelectAndShowLineInfo,line_names_curr,hpts});

end

%% plot anatomy of the jump lines

