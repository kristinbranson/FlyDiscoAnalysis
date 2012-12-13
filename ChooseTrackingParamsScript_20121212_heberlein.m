%% set up path
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';

analysis_protocol = '20121212_non_olympiad_heberlein';
datalocparamsfilestr = 'dataloc_params.txt';
params = {'analysis_protocol',analysis_protocol,...
  'datalocparamsfilestr',datalocparamsfilestr,...
  'settingsdir',settingsdir};
rootdir_fixed = '';
rootoutputdir = '';

expfile = '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20121212_non_olympiad_heberlein/expdirs_ChooseCtraxParameters_20121212.txt';

%% data locations
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

expdirs = importdata(expfile);
expdirs(cellfun(@isempty,expdirs)) = [];

% check that experiments exist
isbad = false(1,numel(expdirs));
for i = 1:numel(expdirs),
  if ~exist(expdirs{i},'dir'),
    isbad(i) = true;
    fprintf('Directory %s does not exist\n',expdirs{i});
  end
end

expdirs(isbad) = [];

%% parameters
nstd_unfixed = 4;
min_counts_unfixed = 100;
dampen_usefixed = false;

%% read metadata

moviefilestr = dataloc_params.moviefilestr;
trxfilestr = dataloc_params.ctraxfilestr;
annfilestr = dataloc_params.annfilestr;

experiments_all = [];
issuccess = false(1,numel(expdirs));
for i = 1:numel(expdirs),
  expdir = expdirs{i};
  res = ReadMetadataFile(fullfile(expdir,dataloc_params.metadatafilestr));
  success = ~metadata.flag_aborted && exist(fullfile(expdir,moviefilestr),'file') && ...
    exist(fullfile(expdir,annfilestr),'file') && ...
    exist(fullfile(expdir,trxfilestr),'file');
  
  issuccess(i) = success;
  if success,
    res.file_system_path = expdir;
    if isempty(experiments_all),
      experiments_all = res;
    else
      experiments_all(end+1) = res; %#ok<SAGROW>
    end
  else
    fprintf('Excluding experiment %s from analysis, missing files or aborted\n',expdir);
  end
end

fprintf('Found %d / %d successful experiments\n',numel(experiments_all),numel(expdirs));
expdirs = expdirs(issuccess);

%% select

weight_order = {'screen_type','screen_reason','rig','bowl','date'};
nexps = 100;

[experiments,nchosen] = ChooseExperimentsCondition(experiments_all,nexps,'weight_order',weight_order);

expdirs = {experiments.file_system_path};
nexpdirs = numel(expdirs);

%% fixed experiments to analyze

if isempty(rootdir_fixed),
  isfixed = false;
else
  expdirs_fixed = dir(fullfile(rootdir_fixed,'*_*'));
  expdirs_fixed = cellfun(@(s)fullfile(rootdir_fixed,s),{expdirs_fixed.name},'UniformOutput',false);
  fixed_trxfilestr = '*fixed*.mat';
  isfixed = ~isempty(expdirs_fixed);
end

%% histogramming parameters

nbins = 100;
prctiles = [0,.001,.0025,.005,.01,.02,.05,.1,.5,1,2,5];
prctiles_bound = [.005,100-.005];
prctiles_bound_jump = [95,100];
prctiles_jump = [95,99,99.5,99.9,99.95,99.99];
prctile_bound_walk = 99.975;

%% load data

%trx = Trx(params{:});

trx = [];
for i = 1:nexpdirs,
  if ~isempty(rootoutputdir),
    [~,experiment_name] = myfileparts(experiments(i).file_system_path);
    expdir = fullfile(rootoutputdir,experiment_name);
  else
    expdir = experiments(i).file_system_path;
  end
  annfile = fullfile(expdir,annfilestr);
  trxcurr = read_ann2(annfile,'trx');
  [trxcurr.expdiri] = deal(i);
  trx = structappend(trx,trxcurr);
end

trx_fixed = [];
trx_unfixed = [];
if isfixed,
  nexpdirs_fixed = numel(expdirs_fixed);
  for i = 1:nexpdirs_fixed,
    tmp = dir(fullfile(expdirs_fixed{i},fixed_trxfilestr));
    datenums = [tmp.datenum];
    [~,order] = sort(datenums);
    trxfile = fullfile(expdirs_fixed{i},tmp(order(end)).name);
    trxcurr = load_tracks(trxfile);
    [trxcurr.expdiri] = deal(i);
    trx_fixed = structappend(trx_fixed,trxcurr);
    trxfile = fullfile(expdirs_fixed{i},trxfilestr);
    trxcurr = load_tracks(trxfile);
    [trxcurr.expdiri] = deal(i);
    trx_unfixed = structappend(trx_unfixed,trxcurr);
  end
end


%% histogram the areas, major, minor axis lengths

a = [trx.a]*2;
b = [trx.b]*2;
area = a.*b*pi;
ecc = b ./ a;
expdiri = nan(size(area));
off = 0;
for i = 1:numel(trx),
  expdiri(off+1:off+numel(trx(i).a)) = trx(i).expdiri;
  off = off + numel(trx(i).a);
end

% choose bins
prctiles_sym = [prctiles,fliplr(100-prctiles)];
lims_a = prctile(a,prctiles_bound);
edges_a = linspace(lims_a(1),lims_a(2),nbins+1);
centers_a = (edges_a(1:end-1)+edges_a(2:end))/2;
lims_b = prctile(b,prctiles_bound);
edges_b = linspace(lims_b(1),lims_b(2),nbins+1);
centers_b = (edges_b(1:end-1)+edges_b(2:end))/2;
lims_area = prctile(area,prctiles_bound);
edges_area = linspace(lims_area(1),lims_area(2),nbins+1);
centers_area = (edges_area(1:end-1)+edges_area(2:end))/2;
lims_ecc = prctile(ecc,prctiles_bound);
edges_ecc = linspace(lims_ecc(1),lims_ecc(2),nbins+1);
centers_ecc = (edges_ecc(1:end-1)+edges_ecc(2:end))/2;

if isfixed,
  
  a_fixed = [trx_fixed.a]*2;
  b_fixed = [trx_fixed.b]*2;
  area_fixed = a_fixed.*b_fixed*pi;
  ecc_fixed = b_fixed ./ a_fixed;
  n_fixed = numel(a_fixed);
  expdiri_fixed = nan(size(area_fixed));
  off = 0;
  for i = 1:numel(trx_fixed),
    expdiri_fixed(off+1:off+numel(trx_fixed(i).a)) = trx_fixed(i).expdiri;
    off = off + numel(trx_fixed(i).a);
  end
  
  a_unfixed = [trx_unfixed.a]*2;
  b_unfixed = [trx_unfixed.b]*2;
  area_unfixed = a_unfixed.*b_unfixed*pi;
  ecc_unfixed = b_unfixed ./ a_unfixed;
  n_unfixed = numel(a_unfixed);
  expdiri_unfixed = nan(size(area_unfixed));
  off = 0;
  for i = 1:numel(trx_unfixed),
    expdiri_unfixed(off+1:off+numel(trx_unfixed(i).a)) = trx_unfixed(i).expdiri;
    off = off + numel(trx_unfixed(i).a);
  end  
  
  
  % what fraction of counts per-bin are removed by fixing? do histogramming
  % relative to sizes for the current movie
  counts_area_unfixed = zeros(1,nbins);
  counts_area_fixed = zeros(1,nbins);
  counts_a_unfixed = zeros(1,nbins);
  counts_a_fixed = zeros(1,nbins);
  counts_b_unfixed = zeros(1,nbins);
  counts_b_fixed = zeros(1,nbins);
  counts_ecc_unfixed = zeros(1,nbins);
  counts_ecc_fixed = zeros(1,nbins);
  
  frac_area_unfixed_perexp = nan(nexpdirs_fixed,nbins);
  frac_area_fixed_perexp = nan(nexpdirs_fixed,nbins);
  frac_a_unfixed_perexp = nan(nexpdirs_fixed,nbins);
  frac_a_fixed_perexp = nan(nexpdirs_fixed,nbins);
  frac_b_unfixed_perexp = nan(nexpdirs_fixed,nbins);
  frac_b_fixed_perexp = nan(nexpdirs_fixed,nbins);
  frac_ecc_unfixed_perexp = nan(nexpdirs_fixed,nbins);
  frac_ecc_fixed_perexp = nan(nexpdirs_fixed,nbins);
  
  for i = 1:nexpdirs_fixed,

    % area
    areacurr = area_unfixed(expdiri_unfixed==i);
    [mu,~,idx] = onedimkmeans(areacurr,2);
    sigma = nan(1,2);
    for j = 1:2,
      sigma(j) = std(areacurr(idx==j),1);
    end
    minv = min(mu'-sigma*nstd_unfixed);
    maxv = max(mu'+sigma*nstd_unfixed);
    [~,centers] = SelectHistEdges(nbins,[minv,maxv],'linear');
    counts_area_unfixed = counts_area_unfixed + hist(areacurr,centers);
    counts = hist(areacurr,centers_area);
    frac_area_unfixed_perexp(i,:) = counts / sum(counts);    
    areacurr = area_fixed(expdiri_fixed==i);
    counts_area_fixed = counts_area_fixed + hist(areacurr,centers);
    counts = hist(areacurr,centers_area);
    frac_area_fixed_perexp(i,:) = counts / sum(counts);
    
    % a
    acurr = a_unfixed(expdiri_unfixed==i);
    [mu,~,idx] = onedimkmeans(acurr,2);
    sigma = nan(1,2);
    for j = 1:2,
      sigma(j) = std(acurr(idx==j),1);
    end
    minv = min(mu'-sigma*nstd_unfixed);
    maxv = max(mu'+sigma*nstd_unfixed);
    [~,centers] = SelectHistEdges(nbins,[minv,maxv],'linear');
    counts_a_unfixed = counts_a_unfixed + hist(acurr,centers);
    counts = hist(acurr,centers_a);
    frac_a_unfixed_perexp(i,:) = counts / sum(counts);    
    acurr = a_fixed(expdiri_fixed==i);
    counts_a_fixed = counts_a_fixed + hist(acurr,centers);
    counts = hist(acurr,centers_a);
    frac_a_fixed_perexp(i,:) = counts / sum(counts);
    
    % b
    bcurr = b_unfixed(expdiri_unfixed==i);
    [mu,~,idx] = onedimkmeans(bcurr,2);
    sigma = nan(1,2);
    for j = 1:2,
      sigma(j) = std(bcurr(idx==j),1);
    end
    minv = min(mu'-sigma*nstd_unfixed);
    maxv = max(mu'+sigma*nstd_unfixed);
    [~,centers] = SelectHistEdges(nbins,[minv,maxv],'linear');
    counts_b_unfixed = counts_b_unfixed + hist(bcurr,centers);
    counts = hist(bcurr,centers_b);
    frac_b_unfixed_perexp(i,:) = counts / sum(counts);    
    bcurr = b_fixed(expdiri_fixed==i);
    counts_b_fixed = counts_b_fixed + hist(bcurr,centers);
    counts = hist(bcurr,centers_b);
    frac_b_fixed_perexp(i,:) = counts / sum(counts);

    % ecc
    ecccurr = ecc_unfixed(expdiri_unfixed==i);
    [mu,~,idx] = onedimkmeans(ecccurr,2);
    sigma = nan(1,2);
    for j = 1:2,
      sigma(j) = std(ecccurr(idx==j),1);
    end
    minv = min(mu'-sigma*nstd_unfixed);
    maxv = max(mu'+sigma*nstd_unfixed);
    [~,centers] = SelectHistEdges(nbins,[minv,maxv],'linear');
    counts_ecc_unfixed = counts_ecc_unfixed + hist(ecccurr,centers);
    counts = hist(ecccurr,centers_ecc);
    frac_ecc_unfixed_perexp(i,:) = counts / sum(counts);    
    ecccurr = ecc_fixed(expdiri_fixed==i);
    counts_ecc_fixed = counts_ecc_fixed + hist(ecccurr,centers);
    counts = hist(ecccurr,centers_ecc);
    frac_ecc_fixed_perexp(i,:) = counts / sum(counts);
    
  end
  
  fix_correction_area = counts_area_fixed ./ counts_area_unfixed;
  fix_correction_area(counts_area_unfixed < min_counts_unfixed) = 1;
  fix_correction_a = counts_a_fixed ./ counts_a_unfixed;
  fix_correction_a(counts_a_unfixed < min_counts_unfixed) = 1;
  fix_correction_b = counts_b_fixed ./ counts_b_unfixed;
  fix_correction_b(counts_b_unfixed < min_counts_unfixed) = 1;
  fix_correction_ecc = counts_ecc_fixed ./ counts_ecc_unfixed;
  fix_correction_ecc(counts_ecc_unfixed < min_counts_unfixed) = 1;
  
  % correct by this fraction
  weight_area = nan(size(area));
  for i = 1:nexpdirs,
    expidx = expdiri==i;
    areacurr = area(expidx);
    [mu,~,idx] = onedimkmeans(areacurr,2);
    sigma = nan(1,2);
    for j = 1:2,
      sigma(j) = std(areacurr(idx==j),1);
    end
    minv = min(mu'-sigma*nstd_unfixed);
    maxv = max(mu'+sigma*nstd_unfixed);
    [~,centers] = SelectHistEdges(nbins,[minv,maxv],'linear');
    [~,~,bin] = myhist(areacurr,centers);
    weight_area(expidx) = fix_correction_area(bin);
  end

  weight_a = nan(size(a));
  for i = 1:nexpdirs,
    expidx = expdiri==i;
    acurr = a(expidx);
    [mu,~,idx] = onedimkmeans(acurr,2);
    sigma = nan(1,2);
    for j = 1:2,
      sigma(j) = std(acurr(idx==j),1);
    end
    minv = min(mu'-sigma*nstd_unfixed);
    maxv = max(mu'+sigma*nstd_unfixed);
    [~,centers] = SelectHistEdges(nbins,[minv,maxv],'linear');
    [~,~,bin] = myhist(acurr,centers);
    weight_a(expidx) = fix_correction_a(bin);
  end

  weight_b = nan(size(b));
  for i = 1:nexpdirs,
    expidx = expdiri==i;
    bcurr = b(expidx);
    [mu,~,idx] = onedimkmeans(bcurr,2);
    sigma = nan(1,2);
    for j = 1:2,
      sigma(j) = std(bcurr(idx==j),1);
    end
    minv = min(mu'-sigma*nstd_unfixed);
    maxv = max(mu'+sigma*nstd_unfixed);
    [~,centers] = SelectHistEdges(nbins,[minv,maxv],'linear');
    [~,~,bin] = myhist(bcurr,centers);
    weight_b(expidx) = fix_correction_b(bin);
  end

  weight_ecc = nan(size(ecc));
  for i = 1:nexpdirs,
    expidx = expdiri==i;
    ecccurr = ecc(expidx);
    [mu,~,idx] = onedimkmeans(ecccurr,2);
    sigma = nan(1,2);
    for j = 1:2,
      sigma(j) = std(ecccurr(idx==j),1);
    end
    minv = min(mu'-sigma*nstd_unfixed);
    maxv = max(mu'+sigma*nstd_unfixed);
    [~,centers] = SelectHistEdges(nbins,[minv,maxv],'linear');
    [~,~,bin] = myhist(ecccurr,centers);
    weight_ecc(expidx) = fix_correction_ecc(bin);
  end

  counts_a_fixed = hist(a_fixed,centers_a);
  %counts_a_fixed = counts_a_fixed(1:end-1);
  frac_a_fixed = counts_a_fixed / n_fixed;
  counts_b_fixed = hist(b_fixed,centers_b);
  %counts_b_fixed = counts_b_fixed(1:end-1);
  frac_b_fixed = counts_b_fixed / n_fixed;
  counts_area_fixed = hist(area_fixed,centers_area);
  %counts_area_fixed = counts_area_fixed(1:end-1);
  frac_area_fixed = counts_area_fixed / n_fixed;
  counts_ecc_fixed = hist(ecc_fixed,centers_ecc);
  %counts_ecc_fixed = counts_ecc_fixed(1:end-1);
  frac_ecc_fixed = counts_ecc_fixed / n_fixed;
  
  counts_a_unfixed = hist(a_unfixed,centers_a);
  %counts_a_unfixed = counts_a_unfixed(1:end-1);
  frac_a_unfixed = counts_a_unfixed / n_unfixed;
  counts_b_unfixed = hist(b_unfixed,centers_b);
  %counts_b_unfixed = counts_b_unfixed(1:end-1);
  frac_b_unfixed = counts_b_unfixed / n_unfixed;
  counts_area_unfixed = hist(area_unfixed,centers_area);
  %counts_area_unfixed = counts_area_unfixed(1:end-1);
  frac_area_unfixed = counts_area_unfixed / n_unfixed;
  counts_ecc_unfixed = hist(ecc_unfixed,centers_ecc);
  %counts_ecc_unfixed = counts_ecc_unfixed(1:end-1);
  frac_ecc_unfixed = counts_ecc_unfixed / n_unfixed;

else
  
  weight_area = ones(size(area));
  weight_a = ones(size(a));
  weight_b = ones(size(b));
  weight_ecc = ones(size(ecc));

end


% histogram data
[counts_a,~,bin_a] = myhist(a,centers_a,'weights',weight_a);
%counts_a = counts_a(1:end-1);
frac_a = counts_a / sum(counts_a);
[counts_b,~,bin_b] = myhist(b,centers_b,'weights',weight_b);
%counts_b = counts_b(1:end-1);
frac_b = counts_b / sum(counts_b);
[counts_area,~,bin_area] = myhist(area,centers_area,'weights',weight_area);
%counts_area = counts_area(1:end-1);
frac_area = counts_area / sum(counts_area);
[counts_ecc,~,bin_ecc] = myhist(ecc,centers_ecc,'weights',weight_ecc);
%counts_ecc = counts_ecc(1:end-1);
frac_ecc = counts_ecc / sum(counts_ecc);

frac_a_perexp = nan(nexpdirs,nbins);
frac_b_perexp = nan(nexpdirs,nbins);
frac_area_perexp = nan(nexpdirs,nbins);
frac_ecc_perexp = nan(nexpdirs,nbins);
mean_area_perexp = nan(1,nexpdirs);

for i = 1:nexpdirs,
  idx = [trx.expdiri]==i;
  
  a_perexp = [trx(idx).a]*2;
  b_perexp = [trx(idx).b]*2;
  area_perexp = a_perexp.*b_perexp*pi;
  ecc_perexp = b_perexp ./ a_perexp;
  n_perexp = numel(a_perexp);
  mean_area_perexp(i) = nansum(area_perexp.*weight_area(expdiri==i))/nansum(weight_area(expdiri==i));
  
  counts_area_perexp = myhist(area_perexp,centers_area,'weights',weight_area(expdiri==i));
  frac_area_perexp(i,:) = counts_area_perexp / sum(counts_area_perexp);
  counts_a_perexp = myhist(a_perexp,centers_a,'weights',weight_a(expdiri==i));
  frac_a_perexp(i,:) = counts_a_perexp / sum(counts_a_perexp);
  counts_b_perexp = myhist(b_perexp,centers_b,'weights',weight_b(expdiri==i));
  frac_b_perexp(i,:) = counts_b_perexp / sum(counts_b_perexp);
  counts_ecc_perexp = myhist(ecc_perexp,centers_ecc,'weights',weight_ecc(expdiri==i));
  frac_ecc_perexp(i,:) = counts_ecc_perexp / sum(counts_ecc_perexp);
  
  tmp = hist(area_perexp,centers_area);
  tmp = tmp / sum(tmp);
  clf;
  maxv = max(max(tmp),max(frac_area_perexp(i,:)));
  plot(centers_area,tmp/maxv,'ko-','markerfacecolor','k');
  hold on;
  plot(centers_area,frac_area_perexp(i,:)/maxv,'rd-','markerfacecolor','r');
  plot(centers_area,frac_area_perexp(i,:)./tmp,'b.-');
  fprintf('experiment %d:\n',i)
  pause(.1);
  
  counts_ecc_perexp = hist(ecc_perexp,centers_ecc,'weights',weight_ecc(expdiri==i));
  %counts_ecc_perexp = counts_ecc_perexp(1:end-1);
  frac_ecc_perexp(i,:) = counts_ecc_perexp / sum(counts_ecc_perexp);
  
end
[~,exporderplot] = sort(abs(mean_area_perexp-median(mean_area_perexp)));

% compute percentiles of data
prctiles_a = weighted_prctile(a,prctiles_sym,weight_a);
prctiles_b = weighted_prctile(b,prctiles_sym,weight_b);
prctiles_area = weighted_prctile(area,prctiles_sym,weight_area);
prctiles_ecc = weighted_prctile(ecc,prctiles_sym,weight_ecc);

%% plot histograms

hfig = 1;
figure(hfig);
clf;

hax = subplot(2,2,1);
colors = jet(nchosen)*.9;
cla;
hold on;
for i = exporderplot,
  plot(centers_area,frac_area_perexp(i,:),'-','color',colors(i,:));
end
if isfixed,
  for i = 1:nexpdirs_fixed,
    hfixed = plot(centers_area,frac_area_fixed_perexp(i,:),'s-','color',[.5,0,0],'markerfacecolor',[.5,0,0]);
    hunfixed = plot(centers_area,frac_area_unfixed_perexp(i,:),'d-','color',[0,.5,.5],'markerfacecolor',[0,.5,.5]);
  end
end
hall = plot(centers_area,frac_area,'ko-','linewidth',3,'markerfacecolor','k');
axisalmosttight;
ylim = get(hax,'YLim');
plot(repmat(prctiles_area,2,1),repmat(ylim',[1,numel(prctiles_sym)]),'g-');
for i = 1:numel(prctiles_sym),
  text(prctiles_area(i),ylim(1)+(ylim(2)-ylim(1))*i/numel(prctiles_sym),[num2str(prctiles_sym(i)),'%'],'color','m');
end
title('area (px^2)');
if isfixed,
  legend([hall,hfixed,hunfixed],{'All exps','Fixed','Pre-fix'},'Location','NorthWest');
else
  legend(hall,{'All exps'},'Location','NorthWest');
end

hax = subplot(2,2,2);
cla;
hold on;
for i = exporderplot,
  plot(centers_a,frac_a_perexp(i,:),'-','color',colors(i,:));
end
if isfixed,
  for i = 1:nexpdirs_fixed,
    hfixed = plot(centers_a,frac_a_fixed_perexp(i,:),'s-','color',[.5,0,0],'markerfacecolor',[.5,0,0]);
    hunfixed = plot(centers_a,frac_a_unfixed_perexp(i,:),'d-','color',[0,.5,.5],'markerfacecolor',[0,.5,.5]);
  end
end
plot(centers_a,frac_a,'ko-','linewidth',3,'markerfacecolor','k');
axisalmosttight;
ylim = get(hax,'YLim');
plot(repmat(prctiles_a,2,1),repmat(ylim',[1,numel(prctiles_sym)]),'g-');
for i = 1:numel(prctiles_sym),
  text(prctiles_a(i),ylim(1)+(ylim(2)-ylim(1))*i/numel(prctiles_sym),[num2str(prctiles_sym(i)),'%'],'color','m');
end
title('semi-major (px)');

hax = subplot(2,2,3);
cla;
hold on;
for i = exporderplot,
  plot(centers_b,frac_b_perexp(i,:),'-','color',colors(i,:));
end
if isfixed,
  for i = 1:nexpdirs_fixed,
    hfixed = plot(centers_b,frac_b_fixed_perexp(i,:),'s-','color',[.5,0,0],'markerfacecolor',[.5,0,0]);
    hunfixed = plot(centers_b,frac_b_unfixed_perexp(i,:),'d-','color',[0,.5,.5],'markerfacecolor',[0,.5,.5]);
  end
end
plot(centers_b,frac_b,'ko-','linewidth',3,'markerfacecolor','k');
axisalmosttight;
ylim = get(hax,'YLim');
plot(repmat(prctiles_b,2,1),repmat(ylim',[1,numel(prctiles_sym)]),'g-');
for i = 1:numel(prctiles_sym),
  text(prctiles_b(i),ylim(1)+(ylim(2)-ylim(1))*i/numel(prctiles_sym),[num2str(prctiles_sym(i)),'%'],'color','m');
end
title('semi-minor (px)');

hax = subplot(2,2,4);
cla;
hold on;
for i = exporderplot,
  plot(centers_ecc,frac_ecc_perexp(i,:),'-','color',colors(i,:));
end
if isfixed,
  for i = 1:nexpdirs_fixed,
    hfixed = plot(centers_ecc,frac_ecc_fixed_perexp(i,:),'s-','color',[.5,0,0],'markerfacecolor',[.5,0,0]);
    hunfixed = plot(centers_ecc,frac_ecc_unfixed_perexp(i,:),'d-','color',[0,.5,.5],'markerfacecolor',[0,.5,.5]);
  end
end
plot(centers_ecc,frac_ecc,'ko-','linewidth',3,'markerfacecolor','k');
axisalmosttight;
ylim = get(hax,'YLim');
plot(repmat(prctiles_ecc,2,1),repmat(ylim',[1,numel(prctiles_sym)]),'g-');
for i = 1:numel(prctiles_sym),
  text(prctiles_ecc(i),ylim(1)+(ylim(2)-ylim(1))*i/numel(prctiles_sym),[num2str(prctiles_sym(i)),'%'],'color','m');
end
title('ecc');

if isfixed,
  
  hfig = 3;
  figure(hfig);
  clf;
  subplot(2,2,1);
  plot(centers_area,fix_correction_area,'k.-');
  axisalmosttight;
  title('Area correction');
  subplot(2,2,2);
  plot(centers_a,fix_correction_a,'k.-');
  axisalmosttight;
  title('Semi-major correction');
  subplot(2,2,3);
  plot(centers_b,fix_correction_b,'k.-');
  axisalmosttight;
  title('Semi-minor correction');
  subplot(2,2,4);
  plot(centers_ecc,fix_correction_ecc,'k.-');
  axisalmosttight;
  title('Ecc correction');
end

%% choose percentiles for size

thresh_prctiles_area = [.02,100-.02];
thresh_prctiles_a = [.1,100-.01];
thresh_prctiles_b = [.01,100-.1];
thresh_prctiles_ecc = [.005,100-.2];

if isfixed,
  min_area = weighted_prctile(area,thresh_prctiles_area(1),weight_area);
  max_area = weighted_prctile(area,thresh_prctiles_area(2),weight_area);
  idx = ~isnan(area);
  mean_area = weighted_mean(area(idx)',weight_area(idx));  
else
  min_area = prctile(area,thresh_prctiles_area(1));
  max_area = prctile(area,thresh_prctiles_area(2));
  mean_area = nanmean(area);
end

if isfixed,
  min_a = weighted_prctile(a,thresh_prctiles_a(1),weight_a);
  max_a = weighted_prctile(a,thresh_prctiles_a(2),weight_a);
  idx = ~isnan(a);
  mean_a = weighted_mean(a(idx)',weight_a(idx));  
else
  min_a = prctile(a,thresh_prctiles_a(1));
  max_a = prctile(a,thresh_prctiles_a(2));
  mean_a = nanmean(a);
end

if isfixed,
  min_b = weighted_prctile(b,thresh_prctiles_b(1),weight_b);
  max_b = weighted_prctile(b,thresh_prctiles_b(2),weight_b);
  idx = ~isnan(b);
  mean_b = weighted_mean(b(idx)',weight_b(idx));  
else
  min_b = prctile(b,thresh_prctiles_b(1));
  max_b = prctile(b,thresh_prctiles_b(2));
  mean_b = nanmean(b);
end

if isfixed,
  min_ecc = weighted_prctile(ecc,thresh_prctiles_ecc(1),weight_ecc);
  max_ecc = weighted_prctile(ecc,thresh_prctiles_ecc(2),weight_ecc);
  idx = ~isnan(ecc);
  mean_ecc = weighted_mean(ecc(idx)',weight_ecc(idx));  
else
  min_ecc = prctile(ecc,thresh_prctiles_ecc(1));
  max_ecc = prctile(ecc,thresh_prctiles_ecc(2));
  mean_ecc = nanmean(ecc);
end

fprintf('min_area = %.3f, mean_area = %.3f, max_area = %.3f\n',min_area,mean_area,max_area);
fprintf('min_a = %.3f, mean_a = %.3f, max_a = %.3f\n',min_a/2,mean_a/2,max_a/2);
fprintf('min_b = %.3f, mean_b = %.3f, max_b = %.3f\n',min_b/2,mean_b/2,max_b/2);
fprintf('min_ecc = %.3f, mean_ecc = %.3f, max_ecc = %.3f\n',min_ecc,mean_ecc,max_ecc);

% 20121212_heberlein
% min_area = 50.027, mean_area = 111.303, max_area = 186.837
% min_a = 3.249, mean_a = 4.998, max_a = 6.731
% min_b = 1.090, mean_b = 1.758, max_b = 2.344
% min_ecc = 0.247, mean_ecc = 0.353, max_ecc = 0.492

% 20120220 housing_CS
% min_area = 66.049, mean_area = 100.472, max_area = 148.277
% min_a = 3.395, mean_a = 4.835, max_a = 6.007
% min_b = 1.328, mean_b = 1.651, max_b = 2.089
% min_ecc = 0.262, mean_ecc = 0.343, max_ecc = 0.497

% 20120110
% min_area = 52.053, mean_area = 82.560, max_area = 133.968
% min_a = 3.476, mean_a = 4.379, max_a = 5.248
% min_b = 1.178, mean_b = 1.496, max_b = 2.790
% min_ecc = 0.250, mean_ecc = 0.344, max_ecc = 0.446

% 20111221
% min_area = 59.601, mean_area = 109.039, max_area = 158.319
% min_a = 2.878, mean_a = 4.972, max_a = 6.160
% min_b = 1.203, mean_b = 1.735, max_b = 2.278
% min_ecc = 0.265, mean_ecc = 0.350, max_ecc = 0.449

% 20110804
% min_area = 67.659, mean_area = 110.744, max_area = 153.596
% min_a = 3.721, mean_a = 4.988, max_a = 6.044
% min_b = 1.327, mean_b = 1.758, max_b = 2.225
% min_ecc = 0.278, mean_ecc = 0.353, max_ecc = 0.409

% 20110725
% min_area = 60.292, mean_area = 109.898, max_area = 155.337
% min_a = 3.631, mean_a = 5.075, max_a = 6.178
% min_b = 1.214, mean_b = 1.715, max_b = 2.120
% min_ecc = 0.259, mean_ecc = 0.339, max_ecc = 0.398

% 20110402
% min_area = 59.261, mean_area = 108.602, max_area = 158.895
% min_a = 3.534, mean_a = 4.964, max_a = 6.260
% min_b = 1.240, mean_b = 1.732, max_b = 2.152
% min_ecc = 0.265, mean_ecc = 0.350, max_ecc = 0.413

% 20110401
% min_area = 58.1, mean_area = 108.1, max_area = 155.9
% min_a = 3.5, mean_a = 5.0, max_a = 6.3
% min_b = 1.2, mean_b = 1.7, max_b = 2.1
% min_ecc = 0.265, mean_ecc = 0.350, max_ecc = 0.413

% 20110326
% min_area = 61.7, mean_area = 107.9, max_area = 151.9
% min_a = 3.5, mean_a = 4.9, max_a = 6.1
% min_b = 1.3, mean_b = 1.7, max_b = 2.2
% min_ecc = 0.274, mean_ecc = 0.350, max_ecc = 0.479

% 20110324
% min_area = 68.5, mean_area = 109.4, max_area = 160.0
% min_a = 3.3, mean_a = 5.0, max_a = 6.3
% min_b = 1.3, mean_b = 1.7, max_b = 2.2
% min_ecc = 0.261, mean_ecc = 0.348, max_ecc = 0.477

% 20110111
%min_area = 56.1, mean_area = 104.7, max_area = 152.7
%min_a = 3.2, mean_a = 4.9, max_a = 6.0
%min_b = 1.2, mean_b = 1.7, max_b = 2.2
%min_ecc = 0.259, mean_ecc = 0.347, max_ecc = 0.476

%% choose the dampening parameters

% min w (xprev_i + w*dx_i - xcurr_i).^2 + (yprev_i + w*dx_i - ycurr_i).^2
% d/dw = 0
% -> 2*(xprev_i + w*dx_i - xcurr_i)*dx + 2*(yprev_i + w*dy_i - ycurr_i)*dy_i = 0
% -> (xprev_i-xcurr_i)*dx_i + (yprev_i-ycurr_i)*dy_i + w*(dx_i^2+dy_i^2) = 0
% -> w = ((xcurr_i-xprev_i)*dx_i + (ycurr_i-yprev_i)*dy_i) / 
%        (dx_i^2 + dy_i^2)
% min w (dtheta_curr_i - w*dtheta_prev_i)^2
% d/dw = 0 
% -> (dtheta_curr_i - w*dtheta_prev_i) * dtheta_prev_i = 0
% -> w = (dtheta_curr_i*dtheta_prev_i)/dtheta_prev_i^2 

if dampen_usefixed,
  trx_dampen = trx_fixed;
else
  trx_dampen = trx;
end

numpos = 0;
denpos = 0;
numangle = 0;
denangle = 0;
for flyidx = 1:numel(trx_dampen),
  x = trx_dampen(flyidx).x;
  y = trx_dampen(flyidx).y;
  theta = trx_dampen(flyidx).theta;
  idx = ~isinf(x) & ~isinf(y) & ~isinf(theta);
  x = x(idx);
  y = y(idx);
  theta = theta(idx);
  dx = diff(x);
  dy = diff(y);
  dtheta = modrange(diff(theta),-pi/2,pi/2);
  numpos = numpos + sum(dx(1:end-1).*dx(2:end)) + sum(dy(1:end-1).*dy(2:end));
  denpos = denpos + sum(dx(1:end-1).^2 + dy(1:end-1).^2);
  numangle = numangle + sum(dtheta(1:end-1).*dtheta(2:end));
  denangle = denangle + sum(dtheta(1:end-1).^2);

end

angle_dampen = numangle ./ denangle;
center_dampen = numpos ./ denpos;

fprintf('Constant velocity center position dampening: %f\n',1-center_dampen);
fprintf('Constant velocity angle dampening: %f\n',1-angle_dampen);

% 20121212_heberlein
% Constant velocity center position dampening: 0.195514
% Constant velocity angle dampening: 0.538359

% 20111221
% Constant velocity center position dampening: 0.168682
% Constant velocity angle dampening: 0.539082

% 20110402
%Constant velocity center position dampening: 0.143002
%Constant velocity angle dampening: 0.498991

% 20110401
%Constant velocity center position dampening: 0.138161
%Constant velocity angle dampening: 0.516704

%Constant velocity center position dampening: 0.142983
%Constant velocity angle dampening: 0.529604

%% choose weight of angle in matching criterion based on error

mean_err_pos = 0;
mean_err_angle = 0;
for flyidx = 1:numel(trx),
  x = trx(flyidx).x;
  y = trx(flyidx).y;
  if numel(x) <= 1,
    continue;
  end
  theta = trx(flyidx).theta;
  idx = ~isinf(x) & ~isinf(y) & ~isinf(theta);
  x = x(idx);
  y = y(idx);
  theta = theta(idx);
  dx = diff(x);
  dy = diff(y);
  dtheta = modrange(diff(theta),-pi/2,pi/2);
  err_pos_curr = (x(1:end-1)+dx*center_dampen - x(2:end)).^2 + ...
    (y(1:end-1)+dy*center_dampen - y(2:end)).^2;
  mean_err_pos = mean_err_pos + sum(err_pos_curr);
  err_angle_curr = modrange( theta(1:end-1)+dtheta*angle_dampen - theta(2:end), -pi/2, pi/2).^2;
  mean_err_angle = mean_err_angle + sum(err_angle_curr);
end

ang_dist_wt = mean_err_pos / mean_err_angle;
fprintf('Weight of angle in matching criterion: %f\n',ang_dist_wt);

% 20121212_heberlein
% Weight of angle in matching criterion: 154.395669

% 20111221
% Weight of angle in matching criterion: 148.266042

% 20110402
% Weight of angle in matching criterion: 133.116340

% 20110401
% Weight of angle in matching criterion: 116.924857

% Weight of angle in matching criterion: 128.595462

%% choose max jump distance


if isfixed,
  trx_jump = trx_fixed;
else
  trx_jump = trx;
end

err = [];
distjump = [];
expdiri_jump = [];
for flyidx = 1:numel(trx_jump),
  fprintf('%d / %d\n',flyidx,numel(trx_jump));
  x = trx_jump(flyidx).x;
  if numel(x) <= 1,
    continue;
  end
  y = trx_jump(flyidx).y;
  theta = trx_jump(flyidx).theta;
  dx = diff(x);
  dy = diff(y);
  dtheta = modrange(diff(theta),-pi/2,pi/2);
  distjump_curr = sqrt(dx.^2+dy.^2);
  distjump = [distjump,distjump_curr]; %#ok<AGROW>
  expdiri_jump = [expdiri_jump,zeros(size(distjump_curr))+trx_jump(flyidx).expdiri]; %#ok<AGROW>

  err_pos_curr = (x(1:end-1)+dx*center_dampen - x(2:end)).^2 + ...
    (y(1:end-1)+dy*center_dampen - y(2:end)).^2;
  err_angle_curr = modrange( theta(1:end-1)+dtheta*angle_dampen - theta(2:end), -pi/2, pi/2).^2;
  err_curr = sqrt(err_pos_curr + ang_dist_wt*err_angle_curr);
  err = [err,err_curr]; %#ok<AGROW>
end

n = numel(err);

lims_err = prctile(err,prctiles_bound_jump);
edges_err = linspace(lims_err(1),lims_err(2),nbins+1);
centers_err = (edges_err(1:end-1)+edges_err(2:end))/2;
counts_err = histc(err,edges_err);
counts_err = [counts_err(1:end-2),counts_err(end-1)+counts_err(end)];
frac_err = counts_err / n;

lims_distjump = prctile(distjump,prctiles_bound_jump);
edges_distjump = linspace(lims_distjump(1),lims_distjump(2),nbins+1);
centers_distjump = (edges_distjump(1:end-1)+edges_distjump(2:end))/2;
counts_distjump = histc(distjump,edges_distjump);
counts_distjump = [counts_distjump(1:end-2),counts_distjump(end-1)+counts_distjump(end)];
frac_distjump = counts_distjump / n;

prctiles_err = prctile(err,prctiles_jump);
prctiles_distjump = prctile(distjump,prctiles_jump);

hfig = 2;
figure(hfig);
clf;

hax = subplot(1,2,1);
plot(centers_err,frac_err,'k.-');
set(hax,'YScale','Log');
ylim = get(hax,'YLim');
hold on;
plot(repmat(prctiles_err,2,1),repmat(ylim',[1,numel(prctiles_jump)]),'g-');
for i = 1:numel(prctiles_jump),
  text(prctiles_err(i),ylim(1)+(ylim(2)-ylim(1))*i/numel(prctiles_jump),[num2str(prctiles_jump(i)),'%'],'color','m');
end
title('err (px/frame)');

hax = subplot(1,2,2);
plot(centers_distjump,frac_distjump,'k.-');
set(hax,'YScale','Log');
ylim = get(hax,'YLim');
hold on;
plot(repmat(prctiles_distjump,2,1),repmat(ylim',[1,numel(prctiles_jump)]),'g-');
for i = 1:numel(prctiles_jump),
  text(prctiles_distjump(i),ylim(1)+(ylim(2)-ylim(1))*i/numel(prctiles_jump),[num2str(prctiles_jump(i)),'%'],'color','m');
end
title('distjump (px/frame)');

%% choose jump thresholds

prctile_thresh_maxjump = 120;
prctile_thresh_minjump = 99.95;

if prctile_thresh_maxjump >= 100,
  max_jump = max(max(err),max(distjump))*prctile_thresh_maxjump/100;
else
  max_jump = prctile([err,distjump],prctile_thresh_maxjump);
end

min_jump = prctile(distjump,prctile_thresh_minjump);
fprintf('Max jump = %f\n',max_jump);
fprintf('Min jump = %f\n',min_jump);

% 20121212_heberlein
% Max jump = 225.031989
% Min jump = 36.841147

% 20110402
% Max jump = 133.212798
% Min jump = 14.222188

% 20110401
% Max jump = 133.212798
% Min jump = 14.222188

% Max jump = 143.152546
% Min jump = 17.006804

%% choose chooseorientations parameters

% cost: 
% dcenter = sqrt(dx^2 + dy^2)
% if dcenter >= min_jump, w = 0; 
% else w = min(max_velocity_angle_weight,velocity_angle_weight*dcenter
% cost = (1-w)*dist(thetacurr,thetaprev) + w*dist(thetacurr,phicurr)

% histogram error in speed as a function 
dcenter = [];
errphi = [];
errprev = [];
expdiri_jump = [];
for flyidx = 1:numel(trx_jump),
  x = trx_jump(flyidx).x;
  y = trx_jump(flyidx).y;
  if numel(x) <= 1,
    continue;
  end
  theta = trx_jump(flyidx).theta;
  dx = diff(x);
  dy = diff(y);
  phi = atan2(dy,dx);
  dcenter_curr = sqrt(dx.^2+dy.^2);
  errphi_curr = abs(modrange(phi-theta(2:end),-pi,pi));
  errprev_curr = abs(modrange(theta(2:end)-theta(1:end-1),-pi,pi));
  dcenter = [dcenter,dcenter_curr]; %#ok<AGROW>
  errphi = [errphi,errphi_curr]; %#ok<AGROW>
  errprev = [errprev,errprev_curr]; %#ok<AGROW>
  expdiri_jump = [expdiri_jump,zeros(size(errphi_curr))+trx_jump(flyidx).expdiri]; %#ok<AGROW>
end

nexpdirs_jump = max(expdiri_jump);

hfig = 4;
figure(hfig);
clf;
idx_notjump = dcenter < min_jump;
[~,centers_dcenter] = SelectHistEdges(nbins,[0,min_jump],'linear');
[~,centers_errphi] = SelectHistEdges(nbins,[0,pi],'linear');
counts = hist3([dcenter(idx_notjump);errphi(idx_notjump)]',{centers_dcenter,centers_errphi});
frac = counts / sum(counts(:));
frac = bsxfun(@rdivide,frac,max(frac,[],1));
imagesc(centers_dcenter([1,nbins]),centers_errphi([1,nbins]),frac);
axis xy;
xlabel('dcenter');
ylabel('err phi');
colormap(logscale_colormap(jet(1000),[.00001,pi]));
fprintf('Click to choose dcenter threshold\n');
[dcenter_thresh,~] = ginput(1);
fprintf('Chose dcenter_thresh = %f\n',dcenter_thresh);

% which frames is each predictor more accurate
idx0 = errphi > errprev;
n1 = nnz(errphi < errprev);
d = min(dcenter,dcenter_thresh);
% score = sum_i [ I(s(i) == 1) log(w*d) + 
%                 I(s(i) == 0) log(1-w*d)
% -> n1 * log(w) + sum_i I(s(i)==0) log(1-w*d)
ntry = 1000;
weights_try = linspace(0,1/dcenter_thresh,ntry);
score = nan(1,ntry);
for i = 1:ntry,
  weight_try = weights_try(i);
  score(i) = n1*log(weight_try) + sum(log(1-weight_try*d(idx0)));
end
[maxscore,i] = max(score);
velocity_angle_weight = weights_try(i);
max_velocity_angle_weight = velocity_angle_weight * dcenter_thresh;

fprintf('velocity_angle_weight = %f\n',velocity_angle_weight);
fprintf('max_velocity_angle_weight = %f\n',max_velocity_angle_weight);

% 20121212_heberlein
% Chose dcenter_thresh = 9.889386
% velocity_angle_weight = 0.043727
% max_velocity_angle_weight = 0.432432

% 20110402
%velocity_angle_weight = 0.030749
%max_velocity_angle_weight = 0.184184

% 20110401
%velocity_angle_weight = 0.030407
%max_velocity_angle_weight = 0.199199