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
nclusters_gt = 10;


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

statfnscurr(statidxremove) = [];
statidxcurr(statidxremove) = [];
datacluster(:,statidxremove) = [];
zdatacluster_norm(:,statidxremove) = [];
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

%% compute pairwise distance between lines, ignoring entries for which either has nan

% L1 distance

lined = zeros(nlinescurr,nlinescurr);
if strcmpi(disttransform,'linearthenlog'),
  zdatacluster_transform = nan(size(zdatacluster_norm));
  idx = abs(zdatacluster_norm)<=linear2loginflectionpt;
  zdatacluster_transform(idx) = zdatacluster_norm(idx);
  zdatacluster_transform(~idx) = (linear2loginflectionpt+...
    log(abs(zdatacluster_norm(~idx))-linear2loginflectionpt+1)).*...
    sign(zdatacluster_norm(~idx));
elseif strcmpi(disttransform,'log'),
  zdatacluster_transform = log(abs(zdatacluster_norm)).*sign(zdatacluster_norm);
else
  zdatacluster_transform = zdatacluster_norm;
end

for linei = 1:nlinescurr,
  for linej = linei+1:nlinescurr,
    dcurr = abs(zdatacluster_transform(linei,:)-zdatacluster_transform(linej,:));
    lined(linei,linej) = nanmean(dcurr);
    lined(linej,linei) = lined(linei,linej);
  end
end
linedvec = squareform(lined,'tovector');

%% cluster lines into nclusters_gt groups using kmeans for groundtruthing

% [idx,c,sumd,d] = mykmeans(zdatacluster_transform,nclusters_gt,'Start','furthestfirst',...
%   'Replicates',100,'EmptyAction','singleton','Display','final');
% 
% % choose the line closest to the center for each cluster
% linesgt_cluster = cell(1,nclusters_gt);
% mind = nan(1,nclusters_gt);
% for i = 1:nclusters_gt,
%   
%   idxcurr = find(idx == i);
%   [mind(i),j] = min(d(idxcurr,i),[],1);
%   j = idxcurr(j);
%   linesgt_cluster{i} = line_names{j};
%   
% end
% 
% [~,lineisgt] = ismember(linesgt_cluster,line_names);
% m = zdatacluster(lineisgt,:);
% 
% % plot the statistics of each cluster
% hfig = 1;
% figure(hfig);
% clf;
% hold on;
% colors = jet(nclusters_gt)*.8;
% for i = 1:nclusters_gt,
%   %boxplot(zdatacluster(idx==i,:),'colors',colors(i,:),'labels',shortstatnames);
%   plot(1:nstatscurr,c(i,:),'o-','Color',colors(i,:),'MarkerFaceColor',colors(i,:));
% end
% set(gca,'XTick',1:nstatscurr,'XTickLabel',statfnscurr);
% htick = rotateticklabel(gca);
% 
% [~,statorder] = ismember(statfnscurr,statfnscurr_ordered);
% 
% ylim = [min(min(c(:)),min(m(:))),max(max(c(:)),max(m(:)))];
% dy = diff(ylim);
% ylim = ylim + dy*.05*[-1,1];
% 
% hfig = 2;
% figure(hfig);
% set(hfig,'Units','pixels','Position',[10,10,1400,550]);
% clf;
% axes('Position',[.05,.075,.7,.9]);
% bar(c(:,statorder));
% hleg = legend(statfnscurr_ordered,'Location','EastOutside');
% set(hleg,'Interpreter','none','Position',[0.76492857142857 0.0660227272727271 0.2055 0.907954545454546]);
% xlabel('Cluster');
% ylabel('Stds from control');
% set(gca,'XLim',[0,nclusters_gt+1],'YLim',ylim);
% box off;
% SaveFigLotsOfWays(hfig,'GTClusterCenters20130920');
% 
% hfig = 3;
% figure(hfig);
% set(hfig,'Units','pixels','Position',[10,10,1400,550]);
% clf;
% axes('Position',[.05,.075,.7,.9]);
% bar(m(:,statorder));
% hleg = legend(statfnscurr_ordered,'Location','EastOutside');
% set(hleg,'Interpreter','none','Position',[0.76492857142857 0.0660227272727271 0.2055 0.907954545454546]);
% xlabel('Cluster');
% ylabel('Stds from control');
% set(gca,'XLim',[0,nclusters_gt+1],'YLim',ylim);
% set(gca,'XTick',1:nclusters_gt,'XTickLabel',shortlinenames(lineisgt));
% box off;
% SaveFigLotsOfWays(hfig,'GTClusterLinesSelected20130920');
% 
% fprintf('%s\n',linesgt_cluster{:});

%% compute pairwise distance between stats, ignoring entries for which either has nan

% 1- abs(correlation coeff)
statd = zeros(nstatscurr,nstatscurr);
for stati = 1:nstatscurr,
  ignorei = isnan(zdatacluster_transform(:,stati));
  for statj = stati+1:nstatscurr,
    ignorecurr = ignorei | isnan(zdatacluster_transform(:,statj));
    if all(ignorecurr),
      statd(stati,statj) = nan;
    else
      r = corrcoef(zdatacluster_transform(~ignorecurr,stati),zdatacluster_transform(~ignorecurr,statj));
      statd(stati,statj) = 1 - abs(r(1,2));
    end
    statd(statj,stati) = statd(stati,statj);
  end
end
statdvec = squareform(statd,'tovector');

%% my version of a clustergram

cgobj = clustergram(zdatacluster_norm',...
  'RowLabels',shortstatnames,...
  'ColumnLabels',shortlinenames,...
  'Standardize','none',...
  'Cluster','all',...
  'RowPDist',statdvec,...
  'ColumnPDist',linedvec,...
  'Linkage','average',...
  'OptimalLeafOrder',true,...
  'ImputeFun',@ClustergramImputeFun);

set(cgobj,'Colormap',myredbluecmap(256));

global LINENAMESSELECTED;
global CLUSTERGRAMORDER;

%% show p-values

idxbigger = pvalue_bigger <= pvalue_smaller;
v = log10(min(2*pvalue_smaller,1));
v(idxbigger) = -log10(min(2*pvalue_bigger(idxbigger),1));

cgobj = clustergram(v(:,statidxcurr)',...
  'RowLabels',shortstatnames,...
  'ColumnLabels',shortlinenames,...
  'Standardize','none',...
  'Cluster','all',...
  'RowPDist',statdvec,...
  'ColumnPDist',linedvec,...
  'Linkage','average',...
  'OptimalLeafOrder',true,...
  'ImputeFun',@ClustergramImputeFun);

set(cgobj,'Colormap',myredbluecmap(256));

% i saved these by dbstopping in clustergram line 1353
global LINEORDER STATORDER;

% show the anatomy data
compartments = fieldnames(linestats.int_manual);
ncompartments = numel(compartments);
anatdata = nan(nlines,ncompartments);
for i = 1:ncompartments,
  anatdata(:,i) = linestats.int_manual.(compartments{i});
end

isexpr = anatdata > 1;
for i = 1:ncompartments,
  ba(:,i) = nanmean(zdatacluster_transform(isexpr(:,i),:),1);
end

dist = pdist(ba','cityblock');
z = linkage(dist,'average');
anatorder = optimalleaforder(z, dist);

%%

hfig = 585;
figure(hfig);
clf;
set(hfig,'Units','pixels','Position',[10 10 2338 840]);
tmp = anatdata;
tmp(isnan(tmp)) = -1;
imagesc(tmp(:,anatorder)');
set(gca,'YTick',1:ncompartments,'YTickLabel',compartments(anatorder))
box off;
w = [0,.25,.55,.75,.9,1]';
colorgreen = [0,.7,0];
coloranat = [
  .75,.75,.75
  bsxfun(@times,w,colorgreen)+bsxfun(@times,1-w,ones(1,3))
  ];
colormap(coloranat);
axis xy;

hfig = 586;
figure(hfig);
clf;
set(hfig,'Units','pixels','Position',[10,10,390 839]);
dendrogram(z,0,'Orientation','left','r',anatorder,'labels',compartments)


%% select and display jumping flies

linenames_jump_avoid0 = {
  'GMR_42E06_AE_01'
  'GMR_68C07_AE_01'
  'GMR_73E12_AE_01'
  'GMR_77E01_AE_01'
  'GMR_82D11_AE_01'
  'GMR_92B02_AE_01'
  };


isin = ismember(line_names,linenames_jump_avoid0);
isout = ~isin;
[linenames_jump_avoid,newisin,coeffs,z] = CompareClusters(zdatacluster_norm,isin,isout,line_names,shortstatnames,...
  'statsplot',{'fractime_jump','dnose2ell'});

linenames_jump_avoid = {
  'GMR_42E06_AE_01'
  'GMR_68C07_AE_01'
  'GMR_73E12_AE_01'
  'GMR_77E01_AE_01'
  'GMR_82D11_AE_01'
  'GMR_92B02_AE_01'
  };

%% select and display clumping flies

linenamesin_clump = {
  'GMR_42F02_AE_01'
  'GMR_26H07_AE_01'
  'GMR_72E09_AE_01'
  'GMR_10A11_AE_01'
  'GMR_70B07_AE_01'
  'GMR_69A10_AE_01'
  'GMR_79F06_AE_01'
  'GMR_61C12_AD_01'
  'GMR_49D01_AE_01'
  'GMR_49E01_AE_01'
  'GMR_12F05_AE_01'
  'GMR_11A11_AE_01'
  'GMR_19H10_AE_01'
  'GMR_10B06_AE_01'
  'GMR_55F12_AE_01'
  'GMR_65A11_AE_01'
  'GMR_43E02_AE_01'
  'GMR_26B01_AE_01'
  'GMR_80G12_AE_01'
  'GMR_68B12_AE_01'
  'GMR_11H05_AE_01'
  'GMR_16D06_AE_01'
  'GMR_32G06_AE_01'
  'GMR_18A01_AE_01'
  'GMR_21B01_AE_01'
  'GMR_25H03_AE_01'
  };

linenamesout_clump = {
  'GMR_65C08_AE_01'
  'GMR_68C09_AE_01'
  'GMR_25G04_AE_01'
  'GMR_40H02_AE_01'
  'GMR_56A11_AE_01'
  'GMR_69F08_AE_01'
  'GMR_12F09_AE_01'
  };

isin = ismember(line_names,linenamesin_clump);

% make sure the flies walk a little bit
stati_walk = find(strcmp(statfnscurr,'fractime_flyany_framewalk'));
controlwalk = mean(setstats.means.fractime_flyany_framewalk(setiscontrol));
isin = isin & datacluster(:,stati_walk)' + controlwalk >= .03;
%isout = ismember(line_names,linenamesout_clump);
isout = ~isin;
[linenames_clump_pivot,newisin,coeffs,z] = CompareClusters(zdatacluster_norm,isin,isout,line_names,shortstatnames);

linenames_clump_pivot = {
  'GMR_79F06_AE_01'
  'GMR_43E02_AE_01'
  'GMR_18A01_AE_01'
  'GMR_49D01_AE_01'
  'GMR_61C12_AD_01'
  'GMR_32G06_AE_01'
  'GMR_26B01_AE_01'
  'GMR_26H07_AE_01'
  'GMR_12F05_AE_01'
  'GMR_11H05_AE_01'
  'GMR_66A06_AE_01'
  'GMR_42F02_AE_01'
  'GMR_82H12_AE_01'
  'GMR_70B07_AE_01'
  'GMR_16D06_AE_01'
  'GMR_19H10_AE_01'
  'GMR_72E09_AE_01'
  'GMR_68B12_AE_01'
  'GMR_21D06_AE_01'
  'GMR_11A11_AE_01'
  'GMR_65A11_AE_01'
  'GMR_21B01_AE_01'
  'GMR_10A11_AE_01'
  'GMR_69A10_AE_01'
  'GMR_80G12_AE_01'
  'GMR_25H03_AE_01'
  'GMR_49E01_AE_01'
  'GMR_10B06_AE_01'
  };

% find some lines that are close but don't pivot
stati_nfliesclose = find(strcmp(statfnscurr,'nflies_close_flyany_frameany'));
stati_velmag = find(strcmp(statfnscurr,'velmag_ctr_flyany_frameany'));
stati_pivot = find(strcmp(statfnscurr,'fractime_flyany_framepivotcenter'));
stati_stop = find(strcmp(statfnscurr,'fractime_flyany_framestop'));

idx = find(zdatacluster_norm(:,stati_nfliesclose) >= min(zdatacluster_norm(isin,stati_nfliesclose)) & ...
  zdatacluster_norm(:,stati_velmag) <= max(zdatacluster_norm(isin,stati_velmag)) & ...
  zdatacluster_norm(:,stati_pivot) <= 0 & ...
  datacluster(:,stati_walk) + controlwalk >= .03);

linenames_clump_nopivot = line_names(idx);

% manually chosen subset

linenames_clump_nopivot = {
  'GMR_26H02_AE_01'
  'GMR_27E09_AE_01'
  'GMR_30B09_AE_01'
  'GMR_41B12_AE_01'
  'GMR_65C08_AE_01'
  'GMR_66B05_AE_01'
  'GMR_91A12_AE_01'
  };

% find some lines that stop a lot but aren't close
idx = find(zdatacluster_norm(:,stati_nfliesclose) <= 0 & ...
  zdatacluster_norm(:,stati_stop) >= min(zdatacluster_norm(isin,stati_stop)) & ...
  datacluster(:,stati_walk) + controlwalk >= .03);

linenames_noclump = line_names(idx);

%%

[meanim,minim,pvalue,maxpvaluesig,nlinesread] = ComputeAverageAnatomyAndPValue(LINENAMESSELECTED);

hfig = 1112;
figure(hfig);
clf;
hax = createsubplots(2,2,.05);
%v = alllines_meanim.*(1-alllines_meanim);
%v(~mask) = inf;

imagesc(max(meanim,[],3)','Parent',hax(1)); axis(hax(1),'image'); colorbar('peer',hax(1));
imagesc(max((alllines_meanim),[],3)','Parent',hax(2)); axis(hax(2),'image'); colorbar('peer',hax(2));
imagesc(max((meanim-alllines_meanim),[],3)','Parent',hax(3)); axis(hax(3),'image'); colorbar('peer',hax(3));
imagesc(min(pvalue,[],3)','Parent',hax(4)); axis(hax(4),'image'); colorbar('peer',hax(4));
colormap(kjetsmooth(256))
set(hax(4),'CLim',[0,.05]);
impixelinfo;
linkaxes(hax);


hfig = 1114;
figure(hfig);
clf;
hax = createsubplots(3,1,.05);

imagesc(max(meanim_noclump,[],3)','Parent',hax(1)); axis(hax(1),'image'); colorbar('peer',hax(1));
imagesc(max(meanim_noclump-alllines_meanim,[],3)','Parent',hax(2)); axis(hax(2),'image'); colorbar('peer',hax(2));
imagesc(min(pvalue_noclump,[],3)','Parent',hax(3)); axis(hax(3),'image'); colorbar('peer',hax(3));
colormap(flipud(kjetsmooth(256)))
set(hax(3),'CLim',[0,.05]);
impixelinfo;
linkaxes(hax);

hfig = 1113;
figure(hfig);
clf;
hax = createsubplots(2,1,.05);

imagesc(max(meanim,[],3)','Parent',hax(1)); axis(hax(1),'image'); colorbar('peer',hax(1));
imagesc(min(pvalue,[],3)','Parent',hax(2)); axis(hax(2),'image'); colorbar('peer',hax(2));
colormap(flipud(kjetsmooth(256)))
set(hax(2),'CLim',[0,.05]);
impixelinfo;
linkaxes(hax);
%%
hfig = 845;
figure(hfig);
clf;


nlinescurr = numel(linenames_clump_pivot);
nc = ceil(sqrt(nlinescurr));
nr = ceil((nlinescurr)/nc);
hax = createsubplots(nr,nc,[.01,.01],hfig);
if numel(hax) > nlinescurr,
  delete(hax(nlinescurr+1:end));
  hax = hax(1:nlinescurr);
end

masklt = pvalue<=.05;
for i = 1:numel(linenames_clump_pivot),
  
  line_name = linenames_clump_pivot{i};
  filename = fullfile(anatomydatadir,sprintf('meanim_%s.mat',line_name));
  if ~exist(filename,'file'),
    
    idxtmp = find(strcmp({imdata.line},line_name));
    qi = [imdata(idxtmp).qi];
    [minqi,j] = min(qi);
    j = idxtmp(j);
    fprintf('%s: using stack %s, qi = %f\n',line_name,imdata(j).name,minqi);
    md = struct;
    md.meanim = loadRaw2StackGreen(imdata(j).raw_file_system_path);
    %md.meanim = md.meanim(:,:,:,2);
    
  else
    md = load(filename,'meanim');    
  end

  md.meanim(~masklt) = 0;
  
  imagesc(max(md.meanim,[],3)','Parent',hax(i));
  axis(hax(i),'image');
  
  title(hax(i),linenames_clump_pivot{i},'Interpreter','none');
  set(hax(i),'XTick',[],'YTick',[]);
  
  colormap(kjetsmooth(256));

  drawnow;
  
end

set(hax,'CLim',[0,1]);
linkaxes(hax);

%% choose all lines that chase more

minnstd = 3;
stati = find(strcmp(statfnscurr,'fractime_flyany_framechase'));
lineis = find(zdatacluster(:,stati)>=minnstd);
linedvec_curr = squareform(lined(lineis,lineis),'tovector');

cgobj2 = clustergram(zdatacluster_norm(lineis,:)',...
  'RowLabels',shortstatnames,...
  'ColumnLabels',shortlinenames(lineis),...
  'Standardize','none',...
  'Cluster','all',...
  'RowPDist',statdvec,...
  'ColumnPDist',linedvec_curr,...
  'Linkage','average',...
  'OptimalLeafOrder',true,...
  'ImputeFun',@ClustergramImputeFun);

set(cgobj2,'Colormap',myredbluecmap(256));

%% female chase hits

line_names_femalechase = {
  'GMR_26F09_AE_01'
  'GMR_26E01_AE_01'
  'GMR_21A01_AE_01'
  'GMR_44D11_AE_01'
  'GMR_72C11_AE_01'
  'GMR_20C08_AE_01'
  'GMR_45F11_AE_01'
  'GMR_51B06_AE_01'
  'GMR_71A09_AE_01'
  'GMR_30G01_AE_01'
  'GMR_65G11_AE_01'
  'GMR_48F12_AE_01'
  'GMR_23A07_AE_01'
  'GMR_24D02_AE_01'
  'GMR_45G01_AE_01'
  'GMR_35C10_AE_01'
  'GMR_35C07_AE_01'
  'GMR_31E02_AE_01'
  'GMR_16F02_AE_01'
  };

ShowLineInfo(line_names_femalechase);
[meanim,minim,pvalue,maxpvaluesig,nlinesread] = ComputeAverageAnatomyAndPValue(line_names_femalechase);

alllines_meanim = load('alllines_anatomy_20130617/meanim.mat');
alllines_meanim = alllines_meanim.meanim;
 
hfig = 1112;
figure(hfig);
clf;
hax = createsubplots(2,2,.02);
%v = alllines_meanim.*(1-alllines_meanim);
%v(~mask) = inf;

imagesc(max(meanim,[],3)','Parent',hax(1)); axis(hax(1),'image'); colorbar('peer',hax(1));
imagesc(max((alllines_meanim),[],3)','Parent',hax(2)); axis(hax(2),'image'); colorbar('peer',hax(2));
imagesc(max((meanim-alllines_meanim),[],3)','Parent',hax(3)); axis(hax(3),'image'); colorbar('peer',hax(3));
imagesc(min(pvalue,[],3)','Parent',hax(4)); axis(hax(4),'image'); colorbar('peer',hax(4));
colormap(kjetsmooth(256))
set(hax(4),'CLim',[0,.05]);
impixelinfo;
linkaxes(hax);
axis(hax,'off');
