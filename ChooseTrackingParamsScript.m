%% set up path
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/misc;
addpath /groups/branson/home/bransonk/tracking/code/JCtrax/filehandling;


%% experiments to analyze

moviefilestr = 'movie.ufmf';
trxfilestr = 'ctrax_results.mat';

experiment_params = struct;
% type of data to analyze
experiment_params.protocol = 'CtraxTest20110202';
% what dates should we analyze
experiment_params.daterange = cell(1,2);
% what lines
experiment_params.linename = '';
% whether the experiments did not start
experiment_params.notstarted = false;
% what files need to exist
experiment_params.subreadfiles = {trxfilestr,moviefilestr};
params = struct2paramscell(experiment_params);

[expdirs,expdir_reads,expdir_writes,experiments,rootreaddir,rootwritedir] = ...
  getExperimentDirs(params{:});

%% parameters of locations of data

trx_params = struct;
% directory within which experiment movies are contained
trx_params.rootreaddir = rootreaddir;
% directory within which experiment tracking & analysis data are
% contained
trx_params.rootwritedir = rootwritedir;
% derived, per-frame measurement directory, relative to experiment directory
trx_params.perframedir = 'perframe';
% base name of file containing processed trajectories
trx_params.trxfilestr = trxfilestr;

file_params = struct;
file_params.moviefilestr = moviefilestr;

%% histogramming parameters

nbins = 100;
prctiles = [.01,.05,.1,.5,1,2,5];
prctiles_bound = [.01,100-.01];
prctiles_bound_jump = [95,100];
prctiles_jump = [95,99,99.5,99.9,99.95,99.99];

%% load data

params = struct2paramscell(trx_params);
trx = Trx(params{:});

for i = 1:numel(expdirs),
  moviefile = fullfile(expdir_reads{i},file_params.moviefilestr);
  matfile = fullfile(expdir_writes{i},trx_params.trxfilestr);
  if ~exist(matfile,'file'),
    continue;
  end
  [~,nframes,fid,vidinfo] = get_readframe_fcn(moviefile);
  fclose(fid);
  trx.AddExpDir(expdirs{i},vidinfo);
end

%% histogram the areas, major, minor axis lengths

a = cell2mat(trx.a)*2;
b = cell2mat(trx.b)*2;
area = a.*b*pi;
ecc = b ./ a;
n = numel(a);

prctiles_sym = [prctiles,fliplr(100-prctiles)];
lims_a = prctile(a,prctiles_bound);
edges_a = linspace(lims_a(1),lims_a(2),nbins+1);
centers_a = (edges_a(1:end-1)+edges_a(2:end))/2;
counts_a = histc(a,edges_a);
counts_a = counts_a(1:end-1);
frac_a = counts_a / n;
lims_b = prctile(b,prctiles_bound);
edges_b = linspace(lims_b(1),lims_b(2),nbins+1);
centers_b = (edges_b(1:end-1)+edges_b(2:end))/2;
counts_b = histc(b,edges_b);
counts_b = counts_b(1:end-1);
frac_b = counts_b / n;
lims_area = prctile(area,prctiles_bound);
edges_area = linspace(lims_area(1),lims_area(2),nbins+1);
centers_area = (edges_area(1:end-1)+edges_area(2:end))/2;
counts_area = histc(area,edges_area);
counts_area = counts_area(1:end-1);
frac_area = counts_area / n;
lims_ecc = prctile(ecc,prctiles_bound);
edges_ecc = linspace(lims_ecc(1),lims_ecc(2),nbins+1);
centers_ecc = (edges_ecc(1:end-1)+edges_ecc(2:end))/2;
counts_ecc = histc(ecc,edges_ecc);
counts_ecc = counts_ecc(1:end-1);
frac_ecc = counts_ecc / n;



prctiles_a = prctile(a,prctiles_sym);
prctiles_b = prctile(b,prctiles_sym);
prctiles_area = prctile(area,prctiles_sym);
prctiles_ecc = prctile(ecc,prctiles_sym);

hfig = 1;
figure(hfig);
clf;

hax = subplot(2,2,1);
bar(centers_area,frac_area);
ylim = get(hax,'YLim');
hold on;
plot(repmat(prctiles_area,2,1),repmat(ylim',[1,numel(prctiles_sym)]),'g-');
for i = 1:numel(prctiles_sym),
  text(prctiles_area(i),ylim(1)+(ylim(2)-ylim(1))*i/numel(prctiles_sym),[num2str(prctiles_sym(i)),'%'],'color','m');
end
title('area (px^2)');

hax = subplot(2,2,2);
bar(centers_a,frac_a);
ylim = get(hax,'YLim');
hold on;
plot(repmat(prctiles_a,2,1),repmat(ylim',[1,numel(prctiles_sym)]),'g-');
for i = 1:numel(prctiles_sym),
  text(prctiles_a(i),ylim(1)+(ylim(2)-ylim(1))*i/numel(prctiles_sym),[num2str(prctiles_sym(i)),'%'],'color','m');
end

title('semi-major (px)');

hax = subplot(2,2,3);
bar(centers_b,frac_b);
ylim = get(hax,'YLim');
hold on;
plot(repmat(prctiles_b,2,1),repmat(ylim',[1,numel(prctiles_sym)]),'g-');
for i = 1:numel(prctiles_sym),
  text(prctiles_b(i),ylim(1)+(ylim(2)-ylim(1))*i/numel(prctiles_sym),[num2str(prctiles_sym(i)),'%'],'color','m');
end
title('semi-minor (px)');

hax = subplot(2,2,4);
bar(centers_ecc,frac_ecc);
ylim = get(hax,'YLim');
hold on;
plot(repmat(prctiles_ecc,2,1),repmat(ylim',[1,numel(prctiles_sym)]),'g-');
for i = 1:numel(prctiles_sym),
  text(prctiles_ecc(i),ylim(1)+(ylim(2)-ylim(1))*i/numel(prctiles_sym),[num2str(prctiles_sym(i)),'%'],'color','m');
end
title('ecc');

%% choose percentiles for size

thresh_prctiles_area = [.01,100-.01];
thresh_prctiles_a = [.02,100-.01];
thresh_prctiles_b = [.01,100-.01];
thresh_prctiles_ecc = [.01,100-.05];

min_area = prctile(area,thresh_prctiles_area(1));
max_area = prctile(area,thresh_prctiles_area(2));
mean_area = nanmean(area);

min_a = prctile(a,thresh_prctiles_a(1));
max_a = prctile(a,thresh_prctiles_a(2));
mean_a = nanmean(a);

min_b = prctile(b,thresh_prctiles_b(1));
max_b = prctile(b,thresh_prctiles_b(2));
mean_b = nanmean(b);

min_ecc = prctile(ecc,thresh_prctiles_ecc(1));
max_ecc = prctile(ecc,thresh_prctiles_ecc(2));
mean_ecc = nanmean(ecc);

fprintf('min_area = %.1f, mean_area = %.1f, max_area = %.1f\n',min_area,mean_area,max_area);
fprintf('min_a = %.1f, mean_a = %.1f, max_a = %.1f\n',min_a/2,mean_a/2,max_a/2);
fprintf('min_b = %.1f, mean_b = %.1f, max_b = %.1f\n',min_b/2,mean_b/2,max_b/2);
fprintf('min_ecc = %.3f, mean_ecc = %.3f, max_ecc = %.3f\n',min_ecc,mean_ecc,max_ecc);

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
