%% set up path

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
addpath /groups/branson/bransonlab/projects/olympiad/anatomy/fileio;
addpath /groups/branson/bransonlab/projects/olympiad/cross_assay/trunk/matlab/kristin;

datafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/CollectedPrimaryPerFrameStatsAndAnatomy20130928.mat';
% imagerydatafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/ImageryData20130725.mat';
% vncdatafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/VNCAnnotations20130905.csv';

anatomydatadir = '/groups/branson/bransonlab/projects/olympiad/AverageAnatomyData20140204';

%% load in data

load(datafile);

%% parameters

maxfraclinesmissingdata = 1;
doanatomyprocessing = false;
statfnset = 'many';
disttransform = 'linearthenlog';
linear2loginflectionpt = 3;
nclusters_gt = 10;

% groupname = 'headdown1';
% linefilename = 'linelists_stephen/Olympiad Head Down - Sheet1.csv';
% groupname = 'nodders1';
% linefilename = 'linelists_stephen/Olympiad Nodders - Sheet1_ver2.csv';
groupname = 'shakers1';
linefilename = 'linelists_stephen/Olympiad Shakers - Sheet1.csv';

%% choose some statistics

ScriptSetStatsToAnalyze;

%% choose some lines

fid = fopen(linefilename,'r');
line_names_curr = {};
while true,
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end
  if isempty(s),
    continue;
  end
  line_names_curr{end+1} = strtrim(s);
end

line_names_curr = unique(line_names_curr);

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

% %% collect this data
% 
% datacluster = nan(nlinescurr,nstatscurr);
% for ii = 1:nstatscurr,
%   i = statidxcurr(ii);
%   datacluster(:,ii) = linestats.normmeans.(statfns{i})(lineidxcurr);
% end
% 
% ncompartments = numel(compartments);
% anatdata = nan(nlinescurr,ncompartments);
% for i = 1:ncompartments,
%   anatdata(:,i) = linestats.int_manual.(compartments{i})(lineidxcurr);
% end
% 
% %% hand-selected correlation removal
% 
% normalizeby = {
% %  'fractime_flyany_framestop'
% %   'velmag_ctr_flyany_frameany'
% %   'dcenter_flyany_frameany'
% %   'dist2wall_flyany_frameany'
%   };
% 
% setiscontrol = strcmp({setstats.metadata.line_name},main_control_line_name);
% setdata = nan(nnz(setiscontrol),nstatscurr);
% for ii = 1:nstatscurr,
%   i = statidxcurr(ii);
%   setdata(:,ii) = setstats.normmeans.(statfns{i})(setiscontrol);
% end
% 
% if isempty(normalizeby),
%   datacluster_norm = datacluster;
%   setdata_norm = setdata;
% else
%   [~,idx] = ismember(normalizeby,statfnscurr);
%   datacluster_norm = nan(size(datacluster));
%   datacluster_norm(:,idx(1)) = datacluster(:,idx(1));
%   j = idx(1);
%   fprintf('%d: %s, average abs val = %f\n',j,statfnscurr{j},mean(abs(datacluster_norm(:,j))));
% 
%   normcoeffs = cell(1,nstatscurr);
%   for ii = 2:numel(idx),
%     
%     is = idx(1:ii-1);
%     j = idx(ii);
%     [normcoeffs{j},~,datacluster_norm(:,j)] = regress(datacluster(:,j),[datacluster(:,is),ones(nlinescurr,1)]);
%     fprintf('%d: %s, average abs val = %f\n',j,statfnscurr{j},mean(abs(datacluster_norm(:,j))));
%     
%   end
% 
%   X = [datacluster(:,idx),ones(nlinescurr,1)];
%   for j = setdiff(1:nstatscurr,idx),
%     idxgood = ~isnan(datacluster(:,j));
%     [normcoeffs{j},~,datacluster_norm(idxgood,j)] = regress(datacluster(idxgood,j),X(idxgood,:));
%     fprintf('%d: %s, average abs val = %f\n',j,statfnscurr{j},mean(abs(datacluster_norm(idxgood,j))));
%   end
%   
%   % also normalize the control set stats so that we can z-score
% 
%   setdata_norm = nan(size(setdata));
% 
%   % first feature
%   j = idx(1);
%   setdata_norm(:,j) = setdata(:,j);
% 
%   % next features
%   for ii = 2:numel(idx),
%     is = idx(1:ii-1);
%     j = idx(ii);
%     pred = [setdata(:,is),ones(nnz(setiscontrol),1)]*normcoeffs{j};
%     setdata_norm(:,j) = setdata(:,j) - pred;
%   end
% 
%   % rest of features
%   X = [setdata(:,idx),ones(nnz(setiscontrol),1)];
%   for j = setdiff(1:nstatscurr,idx),
%     pred = X*normcoeffs{j};
%     setdata_norm(:,j) = setdata(:,j) - pred;
%   end
% 
% end
%   
% % compute mean and standard deviations
% munorm = nanmean(setdata_norm,1);
% signorm = nanstd(setdata_norm,1,1);
% 
% % z-score the line data
% zdatacluster_norm = bsxfun(@rdivide,bsxfun(@minus,datacluster_norm,munorm),signorm);
% 
% % non-normalized version
% mucontrol = nanmean(setdata,1);
% sigcontrol = nanstd(setdata,1,1);
% zdatacluster = bsxfun(@rdivide,bsxfun(@minus,datacluster,mucontrol),sigcontrol);
% 
% % remove nans
% zdatacluster_norm_nonan = zdatacluster_norm;
% zdatacluster_norm_nonan(isnan(zdatacluster_norm)) = 0;
% 
% % 
% % % fill in nans with zeros
% % zdatacluster_norm_nonan = zdatacluster_norm;
% % zdatacluster_norm_nonan(isnan(zdatacluster_norm)) = 0;
% % 
% % % remove some features so that the X matrix is full-rank
% % statidxremove_rank = ismember(statfnscurr,...
% %   {'max_wing_angle_flyany_framewingflick'
% %   'max_absdwing_angle_flyany_framewingflick'
% %   'duration_flyany_framewingflick'
% %   'dcenter_flyany_framewingflick'
% %   'dell2nose_flyany_framewingflick'
% %   'wing_anglel_flyany_frameany'
% %   'fractime_flyany_framechase_notwingextension'
% %   'wing_angle_imbalance_flyany_frameany'}');
% 
% % % here is how I selected features to make X full rank
% % tmp = zdatacluster_norm_nonan;
% % tmp(:,statidxremove_rank) = [];
% % maxrank = rank(tmp);
% % tmpnames = shortstatnames(~statidxremove_rank);
% % idxcanremove = false(1,size(tmp,2));
% % for i = 1:size(tmp,2),
% %   [~,~,r] = regress(tmp(:,i),tmp(:,[1:i-1,i+1:size(tmp,2)]));
% %   if sum(abs(r)) <= .1,
% %     fprintf('%d: %s, regression residual sum = %f\n',i,tmpnames{i},sum(abs(r)));
% %   end
% % 
% %   if rank(tmp(:,[1:i-1,i+1:size(tmp,2)])) == maxrank,
% %     fprintf('%d: %s\n',i,tmpnames{i});
% %     idxcanremove(i) = true;
% %   end
% % end
% 
% %% remove statistics without enough data
% 
% statidxremove = find(sum(isnan(zdatacluster_norm),1) >= nlinescurr*maxfraclinesmissingdata);
% statfnscurr0 = statfnscurr;
% statidxcurr0 = statidxcurr;
% datacluster0 = datacluster;
% zdatacluster_norm0 = zdatacluster_norm;
% 
% statfnscurr(statidxremove) = [];
% statidxcurr(statidxremove) = [];
% datacluster(:,statidxremove) = [];
% zdatacluster_norm(:,statidxremove) = [];
% signorm(statidxremove) = [];
% nstatscurr = numel(statidxcurr);
% 
% %% create short names for stats and lines for plotting
% 
% shortstatnames = statfnscurr;
% shortstatnames = regexprep(shortstatnames,'_flyany','');
% shortstatnames = regexprep(shortstatnames,'^(.*)_fly(.*)_(.*)','$1_$3_$2');
% shortstatnames = regexprep(shortstatnames,'^fractime_frame','fractime_');
% shortstatnames = regexprep(shortstatnames,'^duration_frame','duration_');
% shortstatnames = regexprep(shortstatnames,'_frameany','');
% shortstatnames = regexprep(shortstatnames,'frame','');
% 
% shortlinenames = line_names_curr;
% shortlinenames = regexprep(shortlinenames,'GMR_','R');
% shortlinenames = regexprep(shortlinenames,'_AE_01','');
% shortlinenames = regexprep(shortlinenames,'_AD_01','D');
% 
% %% compute pairwise distance between lines, ignoring entries for which either has nan
% 
% % L1 distance
% 
% lined = zeros(nlinescurr,nlinescurr);
% if strcmpi(disttransform,'linearthenlog'),
%   zdatacluster_transform = nan(size(zdatacluster_norm));
%   idx = abs(zdatacluster_norm)<=linear2loginflectionpt;
%   zdatacluster_transform(idx) = zdatacluster_norm(idx);
%   zdatacluster_transform(~idx) = (linear2loginflectionpt+...
%     log(abs(zdatacluster_norm(~idx))-linear2loginflectionpt+1)).*...
%     sign(zdatacluster_norm(~idx));
% elseif strcmpi(disttransform,'log'),
%   zdatacluster_transform = log(abs(zdatacluster_norm)).*sign(zdatacluster_norm);
% else
%   zdatacluster_transform = zdatacluster_norm;
% end
% 
% for linei = 1:nlinescurr,
%   for linej = linei+1:nlinescurr,
%     dcurr = abs(zdatacluster_transform(linei,:)-zdatacluster_transform(linej,:));
%     lined(linei,linej) = nanmean(dcurr);
%     lined(linej,linei) = lined(linei,linej);
%   end
% end
% linedvec = squareform(lined,'tovector');
% 
% %% compute pairwise distance between stats, ignoring entries for which either has nan
% 
% % 1- abs(correlation coeff)
% statd = zeros(nstatscurr,nstatscurr);
% for stati = 1:nstatscurr,
%   ignorei = isnan(zdatacluster_transform(:,stati));
%   for statj = stati+1:nstatscurr,
%     ignorecurr = ignorei | isnan(zdatacluster_transform(:,statj));
%     if nnz(~ignorecurr) <= 1,
%       statd(stati,statj) = 1;
%     else
%       r = corrcoef(zdatacluster_transform(~ignorecurr,stati),zdatacluster_transform(~ignorecurr,statj));
%       statd(stati,statj) = 1 - abs(r(1,2));
%     end
%     statd(statj,stati) = statd(stati,statj);
%   end
% end
% statdvec = squareform(statd,'tovector');
% 
% %% my version of a clustergram
% 
% cgobj = clustergram(zdatacluster_norm',...
%   'RowLabels',shortstatnames,...
%   'ColumnLabels',shortlinenames,...
%   'Standardize','none',...
%   'Cluster','all',...
%   'RowPDist',statdvec,...
%   'ColumnPDist',linedvec,...
%   'Linkage','average',...
%   'OptimalLeafOrder',true,...
%   'ImputeFun',@ClustergramImputeFun);
% 
% set(cgobj,'Colormap',myredbluecmap(256));
% set(cgobj,'DisplayRange',10);
% 
% global LINENAMESSELECTED;
% global CLUSTERGRAMORDER;

%% compute average anatomy for all lines

matfile = sprintf('meanim_%s_v73.mat',groupname);
if exist(matfile,'file'),
  load(matfile);
else
  [meanim,minim,pvalue,maxpvaluesig,nlinesread] = ComputeAverageAnatomyAndPValue(line_names_curr);
  save(matfile,'-v7.3','meanim','minim','pvalue','maxpvaluesig','nlinesread');
end

figure;
imagesc(max(meanim,[],3)'); 
axis image;
colormap(kjetsmooth(256));
colorbar;
axis off;
set(gcf,'Position',[10 10 1500 900]);
savefig(sprintf('meanim_%s.png',groupname),gcf,'png');

meanim_write = uint16(min(2^16-1,2^16*meanim(:,:,1)));
imwrite(meanim_write,sprintf('meanim_%s.tiff',groupname));
for i = 2:size(meanim,3),
  meanim_write = uint16(min(2^16-1,2^16*meanim(:,:,i)));
  imwrite(meanim_write,sprintf('meanim_%s.tiff',groupname),'writemode','append');
end

figure;
imagesc(min(pvalue,[],3)'); 
axis image;
colormap(flipud(kjetsmooth(256)));
colorbar;
set(gca,'CLim',[0,.05]);
axis off;
set(gcf,'Position',[10 10 1500 900]);
title(sprintf('max pvalue significant: %f',maxpvaluesig));
savefig(sprintf('pvalue_%s.png',groupname),gcf,'png');

pvalue_write = uint16(min(2^16-1,2^16*(1-pvalue(:,:,1))));
imwrite(meanim_write,sprintf('pvalue_%s.tiff',groupname));
for i = 2:size(meanim,3),
  pvalue_write = uint16(min(2^16-1,2^16*(1-pvalue(:,:,i))));
  imwrite(pvalue_write,sprintf('pvalue_%s.tiff',groupname),'writemode','append');
end


%% cache for viewing

imdatafile = 'ImageryData20130824.mat';
maskfile = 'FullBrainMaskSymmetric.mat';
cachedir = '/groups/branson/bransonlab/projects/olympiad/AnatomyCacheData20140209';
alllinematfile = fullfile(anatomydatadir,'alllines_anatomy_20140209/meanim.mat');

for i = 1:numel(line_names_curr),
  CachePerLineAnatomyImages(line_names_curr{i},'imdatafile',imdatafile,'maskfile',maskfile,'cachedir',cachedir,...
    'anatomydir',anatomydatadir);
end

CacheGroupAnatomyImages(line_names_curr,'groupname',groupname,...
  'alllinematfile',alllinematfile,...
  'cachedir',cachedir,...
  'maskfile',maskfile,...
  'anatomydir',anatomydatadir);

ShowLineAnatomy(line_names_curr,...
  'imdatafile',imdatafile,'maskfile',maskfile,'cachedir',cachedir,...
  'groupname',groupname,'maskfile',maskfile,'anatomydir',anatomydatadir);

savefig(sprintf('alllines_%s.png',groupname),gcf,'png');