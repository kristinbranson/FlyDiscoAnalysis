%% set up path

addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
addpath /groups/branson/bransonlab/projects/olympiad/SAGE/MATLABInterface/Trunk;
addpath /groups/branson/bransonlab/projects/olympiad/anatomy/fileio;

datafile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/CollectedPrimaryPerFrameStatsAndAnatomy20130912.mat';

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


%% choose a behavior statistic

statfn = 'fractime_flymale_framechase';
statdir = 'descend';
%minnstds = 45;

stati = find(strcmp(statfnscurr,statfn));
[sortedv,order] = sort(zdatacluster_norm(:,stati),statdir);
hfig = 2;
figure(hfig);
clf;
plot(sortedv,'k.-');
input('Adjust zoom, hit enter when done');
[tmpx,minnstds] = ginput(1);
linei = find(abs(sortedv)<minnstds,1)-1;

% compute anatomy image
line_names_plot = line_names_curr(order(1:linei));
[meanim,pvalue,maxpvaluesig,nlinesread] = ComputeAverageAnatomyAndPValue(line_names_plot);

hfig = 1;
figure(hfig);
clf;
hax = createsubplots(2,1,.05);
axes(hax(1));
imagesc(min(1-meanim,[],3)');
axis image;
colorbar;
axes(hax(2));
imagesc(min(pvalue,[],3)');
axis image;
set(gca,'CLim',[0,.01]);
colormap(flipud(kjetsmooth(256)));
impixelinfo;
colorbar;

filename = sprintf('%s_%s_minnstds%f_%s.mat',shortstatnames{stati},statdir,minnstds,datestr(now,'yyyymmddTHHMMSS'));
save('-v7.3',filename,'meanim','pvalue','maxpvaluesig','nlinesread','line_names_plot','statfn','statdir','minnstds');