%% set up path

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
addpath /groups/branson/bransonlab/projects/olympiad/anatomy/fileio;

datafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/CollectedPrimaryPerFrameStatsAndAnatomy20130928.mat';
% imagerydatafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/ImageryData20130725.mat';
% vncdatafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/VNCAnnotations20130905.csv';

anatomydatadir = '/nobackup/branson/AverageAnatomyData20130618';

%% parameters

maxfraclinesmissingdata = 1;
statfnset = 'many';
disttransform = 'linearthenlog';
linear2loginflectionpt = 3;

%% choose some statistics

ScriptSetStatsToAnalyze;

%% choose some lines

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

line_names_femalefence = {
  'GMR_40F04_AE_01'
  'GMR_49C02_AE_01'
  };

line_names_curr = [line_names_femalechase;line_names_femalefence]';


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
datacluster_alllines = nan(nlines,nstatscurr);
for ii = 1:nstatscurr,
  i = statidxcurr(ii);
  datacluster_alllines(:,ii) = linestats.normmeans.(statfns{i});
end

%ncompartments = numel(compartments);
%anatdata = nan(nlinescurr,ncompartments);
%for i = 1:ncompartments,
%  anatdata(:,i) = linestats.int_manual.(compartments{i})(lineidxcurr);
%end

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
  datacluster_alllines_norm = datacluster_alllines;
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
zdatacluster_alllines_norm = bsxfun(@rdivide,bsxfun(@minus,datacluster_alllines_norm,munorm),signorm);
zcontrol_setdata_norm = bsxfun(@rdivide,bsxfun(@minus,setdata_norm,munorm),signorm);

% non-normalized version
mucontrol = nanmean(setdata,1);
sigcontrol = nanstd(setdata,1,1);
zdatacluster = bsxfun(@rdivide,bsxfun(@minus,datacluster,mucontrol),sigcontrol);
zdatacluster_alllines = bsxfun(@rdivide,bsxfun(@minus,datacluster_alllines,mucontrol),sigcontrol);

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
  
  % control set data
  zcontrol_setdata_transform = nan(size(zcontrol_setdata_norm));
  idx = abs(zcontrol_setdata_norm)<=linear2loginflectionpt;
  zcontrol_setdata_transform(idx) = zcontrol_setdata_norm(idx);
  zcontrol_setdata_transform(~idx) = (linear2loginflectionpt+...
    log(abs(zcontrol_setdata_norm(~idx))-linear2loginflectionpt+1)).*...
    sign(zcontrol_setdata_norm(~idx));
  
  % all lines
  
  zdatacluster_alllines_transform = nan(size(zdatacluster_alllines_norm));
  idx = abs(zdatacluster_alllines_norm)<=linear2loginflectionpt;
  zdatacluster_alllines_transform(idx) = zdatacluster_alllines_norm(idx);
  zdatacluster_alllines_transform(~idx) = (linear2loginflectionpt+...
    log(abs(zdatacluster_alllines_norm(~idx))-linear2loginflectionpt+1)).*...
    sign(zdatacluster_alllines_norm(~idx));
  
elseif strcmpi(disttransform,'log'),
  zdatacluster_transform = log(abs(zdatacluster_norm)).*sign(zdatacluster_norm);
  zdatacluster_alllines_transform = log(abs(zdatacluster_alllines_norm)).*sign(zdatacluster_alllines_norm);
  zcontrol_setdata_transform = log(abs(zcontrol_setdata_norm)).*sign(zcontrol_setdata_norm);

else
  zdatacluster_transform = zdatacluster_norm;  
  zdatacluster_alllines_transform = zdatacluster_alllines_norm;  
  zcontrol_setdata_transform = zcontrol_setdata_norm;

end

for linei = 1:nlinescurr,
  for linej = linei+1:nlinescurr,
    dcurr = abs(zdatacluster_transform(linei,:)-zdatacluster_transform(linej,:));
    lined(linei,linej) = nanmean(dcurr);
    lined(linej,linei) = lined(linei,linej);
  end
end
linedvec = squareform(lined,'tovector');

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

statidxplot = ismember(statfnscurr,statfns_social);

newstatd = statd(statidxplot,statidxplot);
newstatd(isnan(newstatd)) = 1;
newstatd = squareform(newstatd,'tovector');

cgobj = clustergram(zdatacluster_norm(:,statidxplot)',...
  'RowLabels',shortstatnames(statidxplot),...
  'ColumnLabels',shortlinenames,...
  'Standardize','none',...
  'Cluster','all',...
  'RowPDist',newstatd,...
  'ColumnPDist',linedvec,...
  'Linkage','average',...
  'OptimalLeafOrder',true,...
  'ImputeFun',@ClustergramImputeFun);

set(cgobj,'Colormap',myredbluecmap(256),'DisplayRange',10);

global LINENAMESSELECTED;
global CLUSTERGRAMORDER;

%% make a webpage with info for all lines

outfilename = 'FemaleChase/FemaleChaseInfo.html';

ShowLineInfo(line_names_curr,'outfilename',outfilename,'metadata',metadata);

%% t-sne projection


X = zdatacluster_transform;
X(isnan(X)) = 0;
mu_tsne = mean(X,1);
X = bsxfun(@minus,X,mu_tsne);

stati = find(strcmp(statfnscurr,'frame_flyany_framestop'));
ydata = tsne(X,zdatacluster_transform(:,stati));
proj_tsne = ydata;


%% plot the t-sne projection

colorbyfn = 'fractime_flyfemale_framechase';
stati = find(strcmp(statfnscurr,colorbyfn));

interactive = true;

hfig = 10;
figure(hfig);
clf;
set(hfig,'Units','pixels','Position',[376 1028 382 258]);
  
x = ydata(:,1);
y = ydata(:,2);
  
xlim = [min(x),max(x)];
ylim = [min(y),max(y)];
dx = diff(xlim);
dy = diff(ylim);
xlim = xlim+.05*dx*[-1,1];
ylim = ylim+.05*dy*[-1,1];
hold on;
hax = gca;

%plot(projcontrol_pca(:,1),projcontrol_pca(:,2),'kx');
%hold on;

ncolors = 256;
colors = jet(ncolors)*.9;
minz = min(zdatacluster_transform(:,stati));
maxz = max(zdatacluster_transform(:,stati));
colori = max(1,ceil((zdatacluster_transform(:,stati)-minz)/(maxz-minz)*ncolors));

hpts = nan(1,nlinescurr);
for i = 1:nlinescurr,
  %hpts(i) = plot(x(i),y(i),'ko','MarkerFaceColor',colors(colori(i),:));
  hpts(i) = plot(x(i),y(i),'.','Color',colors(colori(i),:));
end
axis square;
axis([xlim,ylim]);

box off;
set(hax,'UserData',struct('selectedpts',[]),'Color','w');

drawnow;
  
if interactive,
  hdata_b1b2a1 = SetUpButtonDown_ReturnPointIndex(hax,x,y,{@ButtonDownFcn_SelectAndShowLineInfo,line_names_curr,hpts});
else
  SaveFigLotsOfWays(hfig,sprintf('FemaleChase/TSNE20131031'));
end

%% t-sne projection for all lines


X = [zdatacluster_alllines_transform;zcontrol_setdata_transform];
X(isnan(X)) = 0;
mu_tsne = mean(X,1);
X = bsxfun(@minus,X,mu_tsne);

stati = find(strcmp(statfnscurr,'frame_flyany_framestop'));
ydata_alllines = tsne(X,zdatacluster_transform(:,stati));
proj_tsne_alllines = ydata_alllines(1:nlines,:);
projcontrol_tsne_alllines = ydata_alllines(nlines+1:end,:);

Scontrol_tsne_alllines = cov(projcontrol_tsne_alllines,1);
mucontrol_tsne_alllines = mean(projcontrol_tsne_alllines,1);

%% plot the t-sne projection for all lines

colorbyfn = 'fractime_flyfemale_framechase';
stati = find(strcmp(statfnscurr,colorbyfn));

interactive = true;

hfig = 11;
figure(hfig);
clf;
set(hfig,'Units','pixels','Position',[376 1028 382 258]);
  
x = ydata_alllines(:,1) - mucontrol_tsne_alllines(1);
y = ydata_alllines(:,2) - mucontrol_tsne_alllines(2);
xcontrol = projcontrol_tsne_alllines(:,1)-mucontrol_tsne_alllines(1);
ycontrol = projcontrol_tsne_alllines(:,2)-mucontrol_tsne_alllines(2);
  
plot(xcontrol,ycontrol,'x','Color',[.4,.4,.4]);
hold on;
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

ncolors = 256;
colors = jet(ncolors)*.9;
minz = min(zdatacluster_alllines_transform(:,stati));
maxz = max(zdatacluster_alllines_transform(:,stati));
colori = max(1,ceil((zdatacluster_alllines_transform(:,stati)-minz)/(maxz-minz)*ncolors));

hpts = nan(1,nlines);
lineidxcurr = ismember(line_names,line_names_curr);

for i = find(lineidxcurr(:)'),
  hpts(i) = plot(x(i),y(i),'o','Color',colors(colori(i),:),'MarkerFaceColor',colors(colori(i),:));
end
for i = find(~lineidxcurr(:)'),
  hpts(i) = plot(x(i),y(i),'.','Color',colors(colori(i),:),'MarkerFaceColor',colors(colori(i),:));
end
axis square;
axis([xlim,ylim]);

box off;
set(hax,'UserData',struct('selectedpts',[]),'Color','w');

drawnow;
  
if interactive,
  hdata_b1b2a1 = SetUpButtonDown_ReturnPointIndex(hax,x,y,{@ButtonDownFcn_SelectAndShowLineInfo,line_names,hpts});
else
  SaveFigLotsOfWays(hfig,sprintf('FemaleChase/TSNEAllLines20131031'));
end

