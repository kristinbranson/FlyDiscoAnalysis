addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;

anatdatafile = 'DescendingNeurons/131029_Lines_usedforsplit.csv';

%% parameters

maxfraclinesmissingdata = 1;
statfnset = 'many';
disttransform = 'linearthenlog';
linear2loginflectionpt = 3;
stringid = '20131029';

%% choose some statistics

ScriptSetStatsToAnalyze;

%% read in anatomy data

fid = fopen(anatdatafile,'r');

if true,
  
s = fgetl(fid);
ss = regexp(s,'"\s*([^"]*)\s*"','tokens');
dnnames = [ss{:}];

ndns = numel(dnnames);
isdnexpr = false(0,ndns);
dnlines = {};

while true,
  
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end

  ss = regexp(s,',','split');
  ss = regexprep(ss,'^"','');
  ss = regexprep(ss,'"$','');
  ss = regexprep(ss,'^''','');
  ss = regexprep(ss,'''$','');
  ss = strtrim(ss);
  if isempty(ss),
    continue;
  end
  isdata = ~cellfun(@isempty,ss);
  if ~any(isdata),
    continue;
  end
  newdnlines = setdiff(ss(isdata),dnlines);
  dnlines = [dnlines,newdnlines];
  [~,idx] = ismember(ss(isdata),dnlines);
  idxdata = find(isdata);
  if size(isdnexpr,1) < numel(dnlines),
    isdnexpr = [isdnexpr;false(numel(dnlines)-size(isdnexpr,1),ndns)];
  end
  isdnexpr(sub2ind([numel(dnlines),ndns],idx,idxdata)) = true;
  
end
  
else

% not sure what is in the first row
s = fgetl(fid);
ss = regexp(s,',','split');
nhitsdn = str2double(ss(2:end-1));

% column titles
s = fgetl(fid);
ss = regexp(s,'"\s*([^"]*)\s*"','tokens');
dnnames = [ss{:}];
ndns = numel(dnnames);
isdnexpr = false(0,ndns);
dnlines = {};

while true,
  
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end

  ss = regexp(s,',','split');
  ss = regexprep(ss,'^"','');
  ss = regexprep(ss,'"$','');
  ss = regexprep(ss,'^''','');
  ss = regexprep(ss,'''$','');
  ss = strtrim(ss);
  if isempty(ss{1}),
    %warning('No line name found for %s, nlines read = %d',s,numel(dnlines));
    continue;
  end
  
  isone = cellfun(@(x) strcmp(x,'1'),ss(2:ndns+1));
  isemp = cellfun(@isempty,ss(2:ndns+1));
  if ~all(isone|isemp),
    warning('some elements of line %d are not empty or one: %s',numel(dnlines)+1,s);
  end
  
  idx = find(strcmp(ss{1},dnlines));
  if ~isempty(idx),
    if any(isdnexpr(idx,:) ~= isone),
      warning('Line %s reannotated, does not match previous annotation(s), or-ing these annotations',ss{1});
      isdnexpr(idx,:) = isdnexpr(idx,:) | isone;
    end
  else
    dnlines{end+1} = ss{1}; %#ok<SAGROW>
    isdnexpr(end+1,:) = isone; %#ok<SAGROW>
  end
  
  olds = s;
  
end

end

fclose(fid);

% remove all vienna lines
dnlines = upper(dnlines);
idx = cellfun(@(x) strcmp(x(1:2),'VT'),dnlines);
dnlines(idx) = [];
isdnexpr(idx,:) = [];

%% lines with DN expression

idx = any(isdnexpr,2);
dnlines = dnlines(idx);
isdnexpr = isdnexpr(idx,:);

%% which lines have we screened
dn_line_names = cellfun(@(x) ['GMR_',x,'_AE_01'],dnlines,'UniformOutput',false);
didscreen = ismember(dn_line_names,line_names);

fprintf('Lines screened:\n');
fprintf('%s\n',dn_line_names{didscreen});

fprintf('Lines not screened:\n');
fprintf('%s\n',dn_line_names{~didscreen});

dn_lines_screened = intersect(dn_line_names,line_names);

line_names_curr = dn_lines_screened;

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

[~,idx] = ismember(line_names_curr,dn_line_names);
isdnexpr_curr = isdnexpr(idx,:);

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
  
elseif strcmpi(disttransform,'log'),
  zdatacluster_transform = log(abs(zdatacluster_norm)).*sign(zdatacluster_norm);
  zcontrol_setdata_transform = log(abs(zcontrol_setdata_norm)).*sign(zcontrol_setdata_norm);

else
  zdatacluster_transform = zdatacluster_norm;  
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

%% compute pairwise distance between lines according to anatomy

anatlined = zeros(nlinescurr,nlinescurr);
for linei = 1:nlinescurr,
  for linej = linei+1:nlinescurr,
    
    dcurr = nnz(isdnexpr_curr(linei,:) ~= isdnexpr_curr(linej,:));
    anatlined(linei,linej) = dcurr;
    anatlined(linej,linei) = dcurr;
  end
end
anatlinedvec = squareform(anatlined,'tovector');

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

statidxplot = ismember(statfnscurr,statfns_few);

cgobj = clustergram(zdatacluster_norm(:,statidxplot)',...
  'RowLabels',shortstatnames(statidxplot),...
  'ColumnLabels',shortlinenames,...
  'Standardize','none',...
  'Cluster','all',...
  'RowPDist',squareform(statd(statidxplot,statidxplot),'tovector'),...
  'ColumnPDist',linedvec,...
  'Linkage','average',...
  'OptimalLeafOrder',true,...
  'ImputeFun',@ClustergramImputeFun);

set(cgobj,'Colormap',myredbluecmap(256),'DisplayRange',10);

global LINENAMESSELECTED;
global CLUSTERGRAMORDER;

%% show the anatomy for this ordering of the lines

doplot = any(isdnexpr_curr,1);
LINEORDER = CLUSTERGRAMORDER.Top;

hfig = 585;
figure(hfig);
clf;
set(hfig,'Units','pixels','Position',[10 10 2338 840]);
tmp = isdnexpr_curr;
imagesc(isdnexpr_curr(LINEORDER,doplot)');
set(gca,'YTick',1:nnz(doplot),'YTickLabel',dnnames(doplot))
box off;
w = [0,.25,.55,.75,.9,1]';
colorgreen = [0,.7,0];
coloranat = [
  .75,.75,.75
  bsxfun(@times,w,colorgreen)+bsxfun(@times,1-w,ones(1,3))
  ];
colormap(coloranat);
axis xy;
set(gca,'YDir','reverse','XTick',[],'TickLength',.0025+[0,0]);
set(gca,'Position',[0.1300    0.0100    0.8600    0.9800]);

%% show lines for each descending neuron

for dni = 1:ndns,
  
  idx = isdnexpr_curr(:,dni);
  if nnz(idx) <= 1,
    continue;
  end

  cgobj = clustergram(zdatacluster_norm(idx,statidxplot)',...
    'RowLabels',shortstatnames(statidxplot),...
    'ColumnLabels',shortlinenames(idx),...
    'Standardize','none',...
    'Cluster','all',...
    'RowPDist',squareform(statd(statidxplot,statidxplot),'tovector'),...
    'ColumnPDist',squareform(lined(idx,idx),'tovector'),...
    'Linkage','average',...
    'OptimalLeafOrder',true,...
    'ImputeFun',@ClustergramImputeFun);
  
  set(cgobj,'Colormap',myredbluecmap(256),'DisplayRange',10);
  
  [hfig,hax] = cgobj.printToFigure();
  
  set(hfig,'Renderer','painters');
  drawnow;
  name = strrep(dnnames{dni},'#','');
  name = regexprep(name,'[^a-zA-Z0-9_]','_');
  SaveFigLotsOfWays(hfig,sprintf('DescendingNeurons/DNClustergram%02d_%s',dni,name));
  
  delete(cgobj);
  
end

%% t-sne projection


X = [zdatacluster_transform;zcontrol_setdata_transform];
X(isnan(X)) = 0;
mu_tsne = mean(X,1);
X = bsxfun(@minus,X,mu_tsne);

stati = find(strcmp(statfnscurr,'frame_flyany_framestop'));
ydata = tsne(X,zdatacluster_transform(:,stati));
proj_tsne = ydata(1:nlinescurr,:);
projcontrol_tsne = ydata(nlinescurr+1:end,:);

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

interactive = false;

colori = nan(1,nlinescurr);
for i = 1:nlinescurr,
  colori(i) = find(isdnexpr_curr(i,:),1);
end

colors = jet(ndns);
markers = '+o*.sd^v><ph';

hfig = 10;
figure(hfig);
clf;
set(hfig,'Units','pixels','Position',[376 1028 382 258]);
  
x = ydata(:,1) - mucontrol_tsne(1);
y = ydata(:,2) - mucontrol_tsne(2);
xcontrol = projcontrol_tsne(:,1)-mucontrol_tsne(1);
ycontrol = projcontrol_tsne(:,2)-mucontrol_tsne(2);
  
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
  
hdn = nan(1,ndns);

hpts = nan(1,nlinescurr);
for i = 1:nlinescurr,
  for j = find(isdnexpr(i,:)),
    %hpts(i) = plot(x(i),y(i),'ko','MarkerFaceColor',colors(colori(i),:));
    markeri = mod(j,numel(markers))+1;
    hpts(i) = plot(x(i),y(i),markers(markeri),'Color',colors(j,:));
    if isnan(hdn(j)),
        hdn(j) = hpts(i);
    end
  end
end
axis square;
axis([xlim,ylim]);

box off;
set(hax,'UserData',struct('selectedpts',[]),'Color','k');

drawnow;
  
if interactive,
  hdata_b1b2a1 = SetUpButtonDown_ReturnPointIndex(hax,x,y,{@ButtonDownFcn_SelectAndShowDNLineInfo,line_names_curr,hpts,isdnexpr_curr,dnnames});
else
  SaveFigLotsOfWays(hfig,sprintf('DescendingNeurons/TSNE20131030'));
end

%% for each neuron, plot the average line stats

for dni = 1:ndns,
  if nnz(isdnexpr_curr(:,dni)) < 2,
    continue;
  end
  hfig = dni;
  figure(dni);
  clf;
  set(hfig,'Units','pixels','Position',[364   656   820   409]);
  %set(gca,'Position',[0.05 0.01 0.94 0.95]);
  
  
  x = zdatacluster_norm(isdnexpr_curr(:,dni),statidxplot);
  xbar = zdatacluster_alllines_norm(:,statidxplot);
  
  n = nnz(isdnexpr_curr(:,dni));
  nbar = nlines;
  mu = nanmean(x,1);
  sig = nanstd(x,1,1)/sqrt(n);
  mubar = nanmean(xbar,1);
  sigbar = nanstd(xbar,1,1)/sqrt(nbar);
  
  g = bsxfun(@plus,[zeros(size(x,1),1);ones(size(xbar,1),1)],...
    1:2:2*nnz(statidxplot));
  xg = [x;xbar];
  l = [shortstatnames(statidxplot);repmat({''},[1,nnz(statidxplot)])];
  
  boxplot(xg(:),g(:),'colors',[1,.5,.5;0,0,0],'labels',l(:),...
    'labelorientation','inline','symbol','.');

  
%   g = repmat(1:nnz(statidxplot),[size(xbar,1),1]);
%   xg = [xbar];
%   l = shortstatnames(statidxplot);
% 
%   boxplot(xg(:),g(:),'colors',[0,0,0],'labels',l(:),...
%     'labelorientation','inline','symbol','.');
  hold on;
%   plot(bsxfun(@plus,1.25:nnz(statidxplot)+.25,rand(nbar,1)*.25-.125),xbar,'.','Color',[.5,.5,.5]);
%   hold on;
%   plot(repmat(1.25:nnz(statidxplot)+.25,[2,1]),[mubar-sigbar;mubar+sigbar],'-','Color',[0,0,0]);
%   plot(1.25:nnz(statidxplot)+.25,mubar,'o','Color','w','MarkerFaceColor','k');
%   plot(repmat(.75:nnz(statidxplot)-.25,[2,1]),[mu-sig;mu+sig],'-','Color',[.8,0,0]);
%   plot(.75:nnz(statidxplot)-.25,mu,'o','Color','k','MarkerFaceColor',[.8,0,0]);
%   plot(bsxfun(@plus,.75:nnz(statidxplot)-.25,rand(n,1)*.25-.125),x,'.','Color',[1,.5,.5]);
  h = plot(bsxfun(@plus,1:2:2*nnz(statidxplot)-.25,rand(n,1)*.25-.125),x,'.','Color',[.8,0,0]);

  h = findall(gca,'-regexp','Tag','(Upper)|(Lower)|(Box)');
  y = get(h,'YData');
  y = [y{:}];
  ylim = [min(-10,min(y)),max(10,max(y))];
  ylim = ylim+[-1,1]*.1*diff(ylim);
  set(gca,'YLim',ylim);
  box off;

  name = strrep(strtrim(dnnames{dni}),'#','');
  name = regexprep(name,'[^a-zA-Z0-9_]','_');
  drawnow;
  SaveFigLotsOfWays(hfig,sprintf('DescendingNeurons/Stats_%02d_%s',dni,name));
end

%% show all lines at once this way

xbar = zdatacluster_alllines_norm(:,statidxplot);
mubar = nanmedian(xbar,1);
sigbar = nanmedian(abs(bsxfun(@minus,xbar,mubar)),1);
 
X = nan(ndns,nnz(statidxplot));

for dni = 1:ndns,
  if nnz(isdnexpr_curr(:,dni)) < 1,
    continue;
  end

  x = zdatacluster_norm(isdnexpr_curr(:,dni),statidxplot);
  
  n = nnz(isdnexpr_curr(:,dni));
  nbar = nlines;
  mu = nanmedian(x,1);
  sig = nanmedian(abs(bsxfun(@minus,x,mu)),1);
  
  X(dni,:) = (mu-mubar)./(sigbar+sig);
  
end

idx = sum(isdnexpr_curr,1)>1;
hfig = 34;
figure(hfig);
clf;
set(hfig,'Units','pixels','Position',[10,10,1722,1366]);
imagesc(X(idx,:));
set(gca,'Position',[ 0.1318    0.1332    0.7354    0.8414]);
set(gca,'CLim',[-5,5]);
colormap(myredbluecmap(256));
colorbar;
s = cell(1,nnz(idx));
idx1 = find(idx);
for ii = 1:numel(idx1),
  i = idx1(ii);
  s{ii} = sprintf('%s (%d)',dnnames{i},nnz(isdnexpr_curr(:,i)));
end
set(gca,'YTick',1:nnz(idx),'YTickLabel',s);
set(gca,'XTick',1:nnz(statidxplot),'XTickLabel',shortstatnames(statidxplot));
rotateticklabel(gca);
box off;

%% show all lines for all dns
 
idx = [];
idxstarts = nan(1,ndns);
for dni = 1:ndns,
  idxstarts(dni) = numel(idx)+1;
  idx = [idx;find(isdnexpr_curr(:,dni))];
end
ncurr = numel(idx);

X = zdatacluster_norm(idx,statidxplot);
% [~,order] = sort(X,1);
% [~,order] = sort(order,1);
% Xprct = (order-1)./(nlinescurr-1);

xbar = zdatacluster_alllines_norm(:,statidxplot);
mubar = nanmedian(xbar,1);
sigbar = nanmedian(abs(bsxfun(@minus,xbar,mubar)),1);

Xnorm = bsxfun(@rdivide,bsxfun(@minus,X,mubar),sigbar);

hfig = 458;
figure(hfig);
clf;
set(hfig,'Units','pixels','Position',[10,10,1722,1366]);
imagesc(Xnorm);
set(gca,'Position',[ 0.1318    0.1332    0.7354    0.8414]);
set(gca,'CLim',[-5,5]);
colormap(myredbluecmap(256));
colorbar;

hold on;
for dni = 2:ndns,
  plot([.5,nnz(statidxplot)+.5],idxstarts(dni)-.5+[0,0],'k-','LineWidth',1);
end

s = cell(1,ndns);
for i = 1:ndns,
  s{i} = sprintf('%s (%d)',dnnames{i},nnz(isdnexpr_curr(:,i)));
end
ytick = (idxstarts-1+[idxstarts(2:end),ncurr+1])/2;
goodidx = any(isdnexpr_curr,1);
set(gca,'YTick',ytick(goodidx),'YTickLabel',s(goodidx));
set(gca,'XTick',1:nnz(statidxplot),'XTickLabel',shortstatnames(statidxplot));
rotateticklabel(gca);
set(gca,'TickLength',[0,0]);
box off;

%% make a webpage with info for all lines

outlineresultsdir = '/Volumes/flyolympiad/Olympiad_Screen/fly_bowl/LineResults';
outolympiaddir = '/Volumes/flyolympiad/Olympiad_Screen';
resultsdir = '/misc/public/Kristin2Hiro/DescendingNeurons';

fid = fopen(fullfile(resultsdir,sprintf('DNInfo%s.html',stringid)),'w');

fprintf(fid,'<html>\n<title>Descending Neuron Information</title>\n<body>\n');
fprintf(fid,'<head>\n');
fprintf(fid,'<style>\n');
fprintf(fid,'table\n');
fprintf(fid,'{\n');
fprintf(fid,'border-collapse:collapse;\n');
fprintf(fid,'}\n');
fprintf(fid,'table, td, th\n');
fprintf(fid,'{\n');
fprintf(fid,'border:1px solid black;\n');
fprintf(fid,'}\n');
fprintf(fid,'</style>\n');
fprintf(fid,'</head>\n');
fprintf(fid,'<body>\n');

fprintf(fid,'<p><h1>Descending Neuron Information 20131029</h1></p>\n');
fprintf(fid,'<p><h3>Median statistics per descending neuron, normalized by median absolute deviations</h3></p>\n');
fprintf(fid,'<p><a href="Stats%s.png"><img src="Stats%s.png" width="1000"/></a></p>\n',stringid,stringid);
fprintf(fid,'<p><a href="StatsAllLines%s.png"><img src="StatsAllLines%s.png" width="1000"/></a></p>\n',stringid,stringid);

fprintf(fid,'<ul>\n');

for dni = 1:ndns,

  idx = isdnexpr_curr(:,dni);
  if ~any(idx),
    continue;
  end
  
  name = strrep(dnnames{dni},'#','');
  name = regexprep(name,'[^a-zA-Z0-9_]','_');
  extraimname = sprintf('DNClustergram%02d_%s.png',dni,name);
  if exist(fullfile(resultsdir,extraimname),'file'),
    extraimnames = {extraimname};
  else
    extraimnames = {};
  end
  extraimname = sprintf('Stats_%02d_%s.png',dni,name);
  if exist(fullfile(resultsdir,extraimname),'file'),
    extraimnames{end+1} = extraimname;
  end
  outfilestr = sprintf('DNInfo%02d_%s.html',dni,name);
  outfilename = fullfile(resultsdir,outfilestr);
  ShowDNLineInfo(line_names_curr(idx),isdnexpr_curr(idx,:),dnnames,...
    'extraimnames',extraimnames,'extraimheights',[500,500],...
    'outfilename',outfilename,...
    'outolympiaddir',outolympiaddir,...
    'outlineresultsdir',outlineresultsdir);
  
  fprintf(fid,'<li><a href="%s">%s</a></li>\n',outfilestr,dnnames{dni});

end

fprintf(fid,'</ul>\n</body>\n</html>\n');
fclose(fid);