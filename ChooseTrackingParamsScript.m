%% set up path
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;
if ispc,  
  settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
  rootdir = 'E:\Data\FlyBowl\bowl_data';
else
  settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
  rootdir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';  
  rootdir_fixed = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/Fixed_AR/EP00005_rc8';
  %rootdir = '/groups/branson/bransonlab/tracking_data/olympiad/FlyBowl/CtraxTest20110315';
end
analysis_protocol = '20110222';
datalocparamsfilestr = 'dataloc_params.txt';
params = {'analysis_protocol',analysis_protocol,...
  'datalocparamsfilestr',datalocparamsfilestr,...
  'settingsdir',settingsdir};

%% data locations
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

%% experiments to analyze

moviefilestr = 'movie.ufmf';
trxfilestr = 'ctrax_results.mat';

experiment_params = struct;
% check for failures
experiment_params.checkflags = true;
% rearing condition
experiment_params.rearing_protocol = 'rc8';
% experimental protocol
experiment_params.experiment_protocol = 'EP0005.xml';
% root data dir
experiment_params.rootdir = rootdir;
% what dates should we analyze
%experiment_params.daterange = cell(1,2);
% what lines
%experiment_params.linename = '';
ngal4 = 20;
ncontrols = 20;
tmpparams = struct2paramscell(experiment_params);

experiments_all = SAGEListBowlExperiments(tmpparams{:});

[experiments,ngal4_chosen,ncontrols_chosen] = choose_expdirs(experiments_all,ngal4,ncontrols);
expdirs = {experiments.file_system_path};
nexpdirs = numel(expdirs);

%% fixed experiments to analyze

expdirs_fixed = dir(fullfile(rootdir_fixed,'*_*'));
expdirs_fixed = cellfun(@(s)fullfile(rootdir_fixed,s),{expdirs_fixed.name},'UniformOutput',false);
fixed_trxfilestr = '*fixed*.mat';
isfixed = ~isempty(expdirs_fixed);

%% histogramming parameters

nbins = 100;
prctiles = [.01,.02,.05,.1,.5,1,2,5];
prctiles_bound = [.005,100-.005];
prctiles_bound_jump = [95,100];
prctiles_jump = [95,99,99.5,99.9,99.95,99.99];

%% load data

%trx = Trx(params{:});

trx = [];
for i = 1:nexpdirs,
  trxfile = fullfile(expdirs{i},dataloc_params.trxfilestr);
  tmp = load(trxfile,'trx');
  [tmp.trx.expdiri] = deal(i);
  trx = structappend(trx,tmp.trx);
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
    tmp = load(trxfile,'trx');
    [tmp.trx.expdiri] = deal(i);
    trx_fixed = structappend(trx_fixed,tmp.trx);
    trxfile = fullfile(expdirs_fixed{i},dataloc_params.trxfilestr);
    tmp = load(trxfile,'trx');
    [tmp.trx.expdiri] = deal(i);
    trx_unfixed = structappend(trx_unfixed,tmp.trx);
  end
end

%% histogram the areas, major, minor axis lengths

a = [trx.a]*2;
b = [trx.b]*2;
area = a.*b*pi;
ecc = b ./ a;
n = numel(a);
expdiri = [trx.expdiri];

% cluster area to get normalization for areas per experiment
area_mu = nan(2,nexpdirs);
area_sigma = nan(2,nexpdirs);
for i = 1:nexpdirs,
  areacurr = area(expdiri==i);
  [area_mu(:,i),~,idx] = onedimkmeans(areacurr,2);
  for j = 1:2,
    area_sigma(j,i) = std(areacurr(idx==1),1);
  end
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

% histogram data

[counts_a,~,bin_a] = myhist(a,centers_a);
%counts_a = counts_a(1:end-1);
frac_a = counts_a / n;
[counts_b,~,bin_b] = myhist(b,centers_b);
%counts_b = counts_b(1:end-1);
frac_b = counts_b / n;
[counts_area,~,bin_area] = myhist(area,centers_area);
%counts_area = counts_area(1:end-1);
frac_area = counts_area / n;
[counts_ecc,~,bin_ecc] = myhist(ecc,centers_ecc);
%counts_ecc = counts_ecc(1:end-1);
frac_ecc = counts_ecc / n;

expdiri = [trx.expdiri];
frac_a_perexp = nan(nexpdirs,nbins);
frac_b_perexp = nan(nexpdirs,nbins);
frac_area_perexp = nan(nexpdirs,nbins);
frac_ecc_perexp = nan(nexpdirs,nbins);
mean_area_perexp = nan(1,nexpdirs);

for i = 1:nexpdirs,
  idx = expdiri==i;
  
  a_perexp = [trx(idx).a]*2;
  b_perexp = [trx(idx).b]*2;
  area_perexp = a_perexp.*b_perexp*pi;
  ecc_perexp = b_perexp ./ a_perexp;
  n_perexp = numel(a_perexp);
  mean_area_perexp(i) = nanmean(area_perexp);
  
  counts_a_perexp = hist(a_perexp,centers_a);
  %counts_a_perexp = counts_a_perexp(1:end-1);
  frac_a_perexp(i,:) = counts_a_perexp / n_perexp;
  counts_b_perexp = hist(b_perexp,centers_b);
  %counts_b_perexp = counts_b_perexp(1:end-1);
  frac_b_perexp(i,:) = counts_b_perexp / n_perexp;
  counts_area_perexp = hist(area_perexp,centers_area);
  %counts_area_perexp = counts_area_perexp(1:end-1);
  frac_area_perexp(i,:) = counts_area_perexp / n_perexp;
  counts_ecc_perexp = hist(ecc_perexp,centers_ecc);
  %counts_ecc_perexp = counts_ecc_perexp(1:end-1);
  frac_ecc_perexp(i,:) = counts_ecc_perexp / n_perexp;
  
end
[~,exporderplot] = sort(abs(mean_area_perexp-median(mean_area_perexp)));

if isfixed,
  
  a_fixed = [trx_fixed.a]*2;
  b_fixed = [trx_fixed.b]*2;
  area_fixed = a_fixed.*b_fixed*pi;
  ecc_fixed = b_fixed ./ a_fixed;
  n_fixed = numel(a_fixed);
  expdiri_fixed = [trx_fixed.expdiri];
  
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
  
  a_unfixed = [trx_unfixed.a]*2;
  b_unfixed = [trx_unfixed.b]*2;
  area_unfixed = a_unfixed.*b_unfixed*pi;
  ecc_unfixed = b_unfixed ./ a_unfixed;
  n_unfixed = numel(a_unfixed);
  expdiri_unfixed = [trx_unfixed.expdiri];
  
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
  
  % what fraction of counts per-bin are removed by fixing?
  fix_correction_area = min(1,frac_area_fixed./frac_area_unfixed);
  fix_correction_a = min(1,frac_a_fixed./frac_a_unfixed);
  fix_correction_b = min(1,frac_b_fixed./frac_b_unfixed);
  fix_correction_ecc = min(1,frac_ecc_fixed./frac_ecc_unfixed);
  
  
  % correct by this fraction
  weight_area = ones(1,numel(area));
  idx = bin_area>=1 & bin_area <= nbins;
  weight_area(idx) = fix_correction_area(bin_area(idx));
  
  weight_a = ones(1,numel(a));
  idx = bin_a>=1 & bin_a <= nbins;
  weight_a(idx) = fix_correction_a(bin_a(idx));
  
  weight_b = ones(1,numel(b));
  idx = bin_b>=1 & bin_b <= nbins;
  weight_b(idx) = fix_correction_b(bin_b(idx));
  
  weight_ecc = ones(1,numel(ecc));
  idx = bin_ecc>=1 & bin_ecc <= nbins;
  weight_ecc(idx) = fix_correction_ecc(bin_ecc(idx));
  
  % use corrected data
  prctiles_a = weighted_prctile(a,prctiles_sym,weight_a);
  prctiles_b = weighted_prctile(b,prctiles_sym,weight_b);
  prctiles_area = weighted_prctile(area,prctiles_sym,weight_area);
  prctiles_ecc = weighted_prctile(ecc,prctiles_sym,weight_ecc);
  
  % correct histograms
  frac_area = frac_area .* fix_correction_area;
  frac_a = frac_a .* fix_correction_a;
  frac_b = frac_b .* fix_correction_b;
  frac_ecc = frac_ecc .* fix_correction_ecc;
  
  frac_area_perexp = bsxfun(@times,frac_area_perexp,fix_correction_area);
  frac_a_perexp = bsxfun(@times,frac_a_perexp,fix_correction_a);
  frac_b_perexp = bsxfun(@times,frac_b_perexp,fix_correction_b);
  frac_ecc_perexp = bsxfun(@times,frac_ecc_perexp,fix_correction_ecc);
  
else

  prctiles_a = prctile(a,prctiles_sym);
  prctiles_b = prctile(b,prctiles_sym);
  prctiles_area = prctile(area,prctiles_sym);
  prctiles_ecc = prctile(ecc,prctiles_sym);

end


%% plot histograms

hfig = 1;
figure(hfig);
clf;

hax = subplot(2,2,1);
colors = cat(1,jet(ngal4_chosen)*.9,repmat(.5,[ngal4_chosen,3]));
cla;
hold on;
for i = exporderplot,
  plot(centers_area,frac_area_perexp(i,:),'-','color',colors(i,:));
end
if isfixed,
  hfixed = plot(centers_area,frac_area_fixed,'s-','color',[.5,0,0],'linewidth',3,'markerfacecolor',[.5,0,0]);
  hunfixed = plot(centers_area,frac_area_unfixed,'d-','color',[0,.5,.5],'linewidth',3,'markerfacecolor',[0,.5,.5]);
end
hall = plot(centers_area,frac_area,'ko-','linewidth',3,'markerfacecolor','k');
axisalmosttight;
ylim = get(hax,'YLim');
plot(repmat(prctiles_area,2,1),repmat(ylim',[1,numel(prctiles_sym)]),'g-');
for i = 1:numel(prctiles_sym),
  text(prctiles_area(i),ylim(1)+(ylim(2)-ylim(1))*i/numel(prctiles_sym),[num2str(prctiles_sym(i)),'%'],'color','m');
end
title('area (px^2)');
legend([hall,hfixed,hunfixed],{'All exps','Fixed','Pre-fix'},'Location','NorthWest');

hax = subplot(2,2,2);
cla;
hold on;
for i = exporderplot,
  plot(centers_a,frac_a_perexp(i,:),'-','color',colors(i,:));
end
if isfixed,
  plot(centers_a,frac_a_fixed,'s-','color',[.5,0,0],'linewidth',3,'markerfacecolor',[.5,0,0]);
  plot(centers_a,frac_a_unfixed,'d-','color',[0,.5,.5],'linewidth',3,'markerfacecolor',[0,.5,.5]);
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
  plot(centers_b,frac_b_fixed,'s-','color',[.5,0,0],'linewidth',3,'markerfacecolor',[.5,0,0]);
  plot(centers_b,frac_b_unfixed,'d-','color',[0,.5,.5],'linewidth',3,'markerfacecolor',[0,.5,.5]);
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
  plot(centers_ecc,frac_ecc_fixed,'s-','color',[.5,0,0],'linewidth',3,'markerfacecolor',[.5,0,0]);
  plot(centers_ecc,frac_ecc_unfixed,'d-','color',[0,.5,.5],'linewidth',3,'markerfacecolor',[0,.5,.5]);
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
  hax = subplot(2,2,1);
  plot(centers_area,fix_correction_area,'k.-');
  axisalmosttight;
  title('Area correction');
  hax = subplot(2,2,2);
  plot(centers_a,fix_correction_a,'k.-');
  axisalmosttight;
  title('Semi-major correction');
  hax = subplot(2,2,3);
  plot(centers_b,fix_correction_b,'k.-');
  axisalmosttight;
  title('Semi-minor correction');
  hax = subplot(2,2,4);
  plot(centers_ecc,fix_correction_ecc,'k.-');
  axisalmosttight;
  title('Ecc correction');
end

%% choose percentiles for size

thresh_prctiles_area = [.02,100-.01];
thresh_prctiles_a = [.05,100-.01];
thresh_prctiles_b = [.01,100-.05];
thresh_prctiles_ecc = [.01,100-.05];

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

fprintf('min_area = %.1f, mean_area = %.1f, max_area = %.1f\n',min_area,mean_area,max_area);
fprintf('min_a = %.1f, mean_a = %.1f, max_a = %.1f\n',min_a/2,mean_a/2,max_a/2);
fprintf('min_b = %.1f, mean_b = %.1f, max_b = %.1f\n',min_b/2,mean_b/2,max_b/2);
fprintf('min_ecc = %.3f, mean_ecc = %.3f, max_ecc = %.3f\n',min_ecc,mean_ecc,max_ecc);

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
numpos = 0;
denpos = 0;
numangle = 0;
denangle = 0;
for flyidx = 1:trx.nflies,
  x = trx(flyidx).x;
  y = trx(flyidx).y;
  theta = trx(flyidx).theta;
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

%Constant velocity center position dampening: 0.142983
%Constant velocity angle dampening: 0.529604

%% choose weight of angle in matching criterion based on error

mean_err_pos = 0;
mean_err_angle = 0;
for flyidx = 1:trx.nflies,
  x = trx(flyidx).x;
  y = trx(flyidx).y;
  theta = trx(flyidx).theta;
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
% Weight of angle in matching criterion: 128.595462

%% choose max jump distance

err = [];
distjump = [];
for flyidx = 1:trx.nflies,
  x = trx(flyidx).x;
  y = trx(flyidx).y;
  theta = trx(flyidx).theta;
  dx = diff(x);
  dy = diff(y);
  dtheta = modrange(diff(theta),-pi/2,pi/2);
  distjump_curr = sqrt(dx.^2+dy.^2);
  distjump = [distjump,distjump_curr];
  err_pos_curr = (x(1:end-1)+dx*center_dampen - x(2:end)).^2 + ...
    (y(1:end-1)+dy*center_dampen - y(2:end)).^2;
  err_angle_curr = modrange( theta(1:end-1)+dtheta*angle_dampen - theta(2:end), -pi/2, pi/2).^2;
  err_curr = sqrt(err_pos_curr + ang_dist_wt*err_angle_curr);
  err = [err,err_curr];
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
prctile_thresh_minjump = 99.925;

if prctile_thresh_maxjump >= 100,
  max_jump = max(max(err),max(distjump))*prctile_thresh_maxjump/100;
else
  max_jump = prctile([err,distjump],prctile_thresh_maxjump);
end

min_jump = prctile(distjump,prctile_thresh_minjump);
fprintf('Max jump = %f\n',max_jump);
fprintf('Min jump = %f\n',min_jump);

% Max jump = 143.152546
% Min jump = 17.006804

%% choose chooseorientations parameters