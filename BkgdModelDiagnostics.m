function [bkgd_diagnostics,res] = BkgdModelDiagnostics(expdir,varargin)

bkgd_diagnostics = struct;

%% parse parameters
[analysis_protocol,settingsdir,datalocparamsfilestr] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'datalocparamsfilestr','dataloc_params.txt');

%% read parameters

datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
bkgdmodeldiagnosticsparamsfile = fullfile(settingsdir,analysis_protocol,...
  dataloc_params.bkgdmodeldiagnosticsparamsfilestr);
params = ReadParams(bkgdmodeldiagnosticsparamsfile);

%% set up figure

res = struct;
hfig = 1;
figure(hfig);
clf(hfig,'reset');
set(hfig,'Units','Pixels','Position',params.figpos);
nax_r = 4;
nax_c = 4;
hax = createsubplots(nax_r,nax_c,.04,hfig);
hax = reshape(hax(1:nax_r*nax_c),[nax_r,nax_c]);
res.hax = hax;
res.hfig = hfig;
imax = [];

%% read annotation file

annfile = fullfile(expdir,dataloc_params.annfilestr);
ann = read_ann(annfile);
% resize images
annfile_images = {'background_median','background_mean',...
  'background_mad','background_std',...
  'fracframesisback',...
  'background_center','background_dev',...
  'isarena'};
nr = ann.movie_width;
nc = ann.movie_height;

% resize images read from annotation
for i = 1:numel(annfile_images),
  fn = annfile_images{i};
  if isfield(ann,fn),
    ann.(fn) = permute(reshape(ann.(fn),[nc,nr,numel(ann.(fn))/(nr*nc)]),[2,1,3]);
  end
end

%% image of background center
bkgd_diagnostics.background_center = ann.background_center;
res.him_bkgdcenter = image(repmat(uint8(ann.background_center),[1,1,3]),'Parent',hax(1,1));
res.hti_bkgdcenter = title(hax(1,1),'Background Center');
axis(hax(1,1),'image','off');
imax(end+1) = hax(1,1);

%% histogram of background center intensity
[edges_bkgdcenter,centers_bkgdcenter] = ...
  get_bins(params.nbins_bkgdcenter,params.lim_bkgdcenter,false,false);
counts = histc(ann.background_center(:)',edges_bkgdcenter);
counts = [counts(1:end-2),counts(end-1)+counts(end)];
frac = counts / numel(ann.background_center);

res.hline_histbkgdcenter = ...
  semilogy_with_zeros(hax(2,1),centers_bkgdcenter,frac,...
  [edges_bkgdcenter(1),edges_bkgdcenter(end)],1,...
  'k.-',params.linestyleparams_bkgdcenter{:});
res.hti_histbkgdcenter = title(hax(2,1),'Bkgd Ctr Px Intensities');
bkgd_diagnostics.frac_bkgdcenter = frac;
bkgd_diagnostics.edges_bkgdcenter = edges_bkgdcenter;

%% background center intensity summary statistics

% bin with highest count
[~,bin] = max(frac);
bkgd_diagnostics.histmode_bkgdcenter = centers_bkgdcenter(bin);

% mean, std, prctiles
bkgd_diagnostics.mean_bkgdcenter = nanmean(ann.background_center(:));
bkgd_diagnostics.std_bkgdcenter = nanstd(ann.background_center(:),1);
bkgd_diagnostics.prctiles_bkgdcenter = prctile(ann.background_center(:),params.prctiles_compute);

%% image of background deviation
dev = [];
devtype = '';
if isfield(ann,'background_mad'),
  dev = ann.background_mad;
  devtype = 'MAD';
elseif isfield(ann,'background_std'),
  dev = ann.background_std;
  devtype = 'Std';
elseif isfield(ann,'background_dev'),
  dev = ann.background_dev;
  devtype = 'Dev';
end

if ~isempty(dev),
  if isnan(params.clim_bkgddev(2)),
    params.clim_bkgddev(2) = max(dev(:));
  end
  if isnan(params.clim_bkgddev(1)),
    params.clim_bkgddev(1) = 0;
  end
  res.him_bkgddev = imagesc(dev,'Parent',hax(1,2),params.clim_bkgddev);
  res.hti_bkgddev = title(hax(1,2),['Background ',devtype]);
  res.hcb_bkgddev = colorbar('peer',hax(1,2),'Location','East');
  axis(hax(1,2),'image','off');
  imax(end+1) = hax(1,2);
  cbpos = get(res.hcb_bkgddev,'Position');
  axpos = get(hax(1,2),'Position');
  cbpos(1) = axpos(1)+axpos(3)+cbpos(3);
  set(res.hcb_bkgddev,'Position',cbpos);
  
  %% histogram of background deviation
  [edges_bkgddev,centers_bkgddev] = get_bins(params.nbins_bkgddev,params.lim_bkgddev,false,true);
  counts = histc(dev(:)',edges_bkgddev);
  counts = [counts(1:end-2),counts(end-1)+counts(end)];
  frac = counts / numel(dev);
  res.hline_histbkgddev = ...
    semilogy_with_zeros(hax(2,2),centers_bkgddev,frac,...
    [edges_bkgddev(1),edges_bkgddev(end)],1,...
    'k.-',params.linestyleparams_bkgddev{:});
  res.hti_histbkgddev = title(hax(2,2),['Bkgd ',devtype,' Hist']);

  bkgd_diagnostics.frac_bkgddev = frac;
  bkgd_diagnostics.edges_bkgddev = edges_bkgddev;

  %% mean, std, mode, percentiles of background median absolute deviation
  [~,bin] = max(frac);
  bkgd_diagnostics.histmode_bkgddev = centers_bkgddev(bin);
  
  % mean, std, prctiles
  bkgd_diagnostics.mean_bkgddev = nanmean(dev(:));
  bkgd_diagnostics.std_bkgddev = nanstd(dev(:),1);
  bkgd_diagnostics.prctiles_bkgddev = prctile(dev(:),params.prctiles_compute);
  
else
  
  bkgd_diagnostics.frac_bkgddev = [];
  bkgd_diagnostics.edges_bkgddev = [];
  bkgd_diagnostics.mean_bkgddev = nan;
  bkgd_diagnostics.std_bkgddev = nan;
  bkgd_diagnostics.histmode_bkgddev = nan;
  bkgd_diagnostics.prctiles_bkgddev = nan(size(params.prctiles_compute));
  
end

%% Image of background center overlaid with areas that are not always
% background (according to the alwaysbkgd mask) in red.
if isfield(ann,'isarena'),
  isarena = ann.isarena~=0;
else
  isarena = true([nr,nc]);
  if all(isfield(ann,{'max_nonarena','background_center'})),
    isarena = isarena & (ann.background_center > ann.max_nonarena);
  end
  if all(isfield(ann,{'min_nonarena','background_center'})),
    isarena = isarena & (ann.background_center < ann.min_nonarena);
  end
  if all(isfield(ann,{'do_set_circular_arena',...
      'arena_center_x','arena_center_y','arena_radius'})) && ...
      ann.do_set_circular_arena ~= 0,
    [gridx,gridy] = meshgrid(0:nc-1,0:nr-1);
    isarena = isarena & ...
      ( (gridx - ann.arena_center_x).^2 + ...
      (gridy - ann.arena_center_y).^2 ) <= ann.arena_radius^2;
  end
end
tmp = ann.background_center;
tmp(~isarena) = 255;
imalwaysbkgd = uint8(cat(3,tmp,repmat(ann.background_center,[1,1,2])));
res.him_alwaysbkgd = image(imalwaysbkgd,'Parent',hax(1,3));
res.hti_alwaysbkgd = title(hax(1,3),'Always Bkgd');
axis(hax(1,3),'image','off');
if all(isarena(:)),
  res.hte_alwaysbkgd = text(nc/2,nr/2,'(None)','color','r',...
    'Parent',hax(1,3),'HorizontalAlignment','center');
end
imax(end+1) = hax(1,3);

bkgd_diagnostics.isarena = isarena;
bkgd_diagnostics.imalwaysbkgd = imalwaysbkgd;

%% fraction of pixels that are always background
bkgd_diagnostics.fracpx_alwaysbkgd = nnz(isarena) / (nr*nc);

%% histogram of background center intensities for pixels that are always
% background
if any(~isarena(:)),
  counts = histc(ann.background_center(~isarena)',edges_bkgdcenter);
  counts = [counts(1:end-2),counts(end-1)+counts(end)];
  frac = counts / nnz(~isarena);
else
  frac = zeros(1,params.nbins_bkgdcenter);
end
res.hline_histalwaysbkgd = ...
  semilogy_with_zeros(hax(2,3),centers_bkgdcenter,frac,...
  [edges_bkgdcenter(1),edges_bkgdcenter(end)],1,...
  'k.-',params.linestyleparams_histalwaysbkgd{:});
res.hti_histalwaysbkgd = title(hax(2,3),'Always Bkgd Px Intensities');

bkgd_diagnostics.frac_alwaysbkgd = frac;

%% always background intensity summary statistics

% bin with highest count
[~,bin] = max(frac);
bkgd_diagnostics.histmode_alwaysbkgd = centers_bkgdcenter(bin);

% mean, std, prctiles
bkgd_diagnostics.mean_alwaysbkgd = nanmean(ann.background_center(~isarena));
bkgd_diagnostics.std_alwaysbkgd = nanstd(ann.background_center(~isarena),1);
bkgd_diagnostics.prctiles_alwaysbkgd = prctile(ann.background_center(~isarena),params.prctiles_compute);

%% image of fraction of frames that are used to estimate background model
if ~isfield(ann,'use_expbgfgmodel') || ...
    ann.use_expbgfgmodel == 0 || ...
    ~isfield(ann,'fracframesisback'),
  fracframesisback = ones(nr,nc);
else
  fracframesisback = ann.fracframesisback;
end
res.him_fracframesisback = imagesc(fracframesisback,'Parent',hax(3,1),[0,1]);
res.hti_fracframesisback = title(hax(3,1),'Frac. Frames Used in Bkgd Model');
axis(hax(3,1),'image','off');
res.hcb_fracframesisback = colorbar('peer',hax(3,1),'Location','East');
cbpos = get(res.hcb_fracframesisback,'Position');
axpos = get(hax(3,1),'Position');
cbpos(1) = axpos(1)+axpos(3)+cbpos(3);
set(res.hcb_fracframesisback,'Position',cbpos);
imax(end+1) = hax(3,1);

bkgd_diagnostics.fracframesisback = fracframesisback;

%% histogram of fraction of frames used to estimate background
[edges_fracframesisback,centers_fracframesisback] = ...
  get_bins(params.nbins_fracframesisback,params.lim_fracframesisback,false,false);
counts = histc(fracframesisback(:)',edges_fracframesisback);
counts = [counts(1:end-2),counts(end-1)+counts(end)];
frac = counts / numel(fracframesisback);
res.hline_histfracframesisback = ...
  semilogy_with_zeros(hax(4,1),centers_fracframesisback,frac,...
  [edges_fracframesisback(1),edges_fracframesisback(end)],1,...
  'k.-',params.linestyleparams_histfracframesisback{:});
res.hti_histfracframesisback = title(hax(4,1),'Frac. Frames Used in Bkgd Model Hist');

bkgd_diagnostics.frac_fracframesisback = frac;
bkgd_diagnostics.edges_fracframesisback = edges_fracframesisback;

%% mean, std, prctiles, mode fraction of frames used to estimate background model

% bin with highest count
[~,bin] = max(frac);
bkgd_diagnostics.histmode_fracframesisback = centers_fracframesisback(bin);

bkgd_diagnostics.mean_fracframesisback = mean(fracframesisback(:));
bkgd_diagnostics.std_fracframesisback = std(fracframesisback(:),1);
bkgd_diagnostics.prctiles_fracframesisback = prctile(fracframesisback(:),params.prctiles_compute);

%% Image of background center overlaid with areas that were interpolated in red
if isfield(ann,'min_frac_frames_isback'),
  isinterpolated = fracframesisback < ann.min_frac_frames_isback;
else
  isinterpolated = false(nr,nc);
end
b = bwboundaries(imdilate(isinterpolated,ones(2)),8);
res.him_isinterpolated = image(repmat(uint8(ann.background_center),[1,1,3]),'Parent',hax(3,3));
hold(hax(3,3),'on');
res.hli_isinterpolated = nan(size(b));
for i = 1:numel(b),
  res.hli_isinterpolated(i) = plot(hax(3,3),b{i}(:,1),b{i}(:,2),'r-');
end

res.hti_isinterpolated = title(hax(3,3),'Bkgd Interpolated');
axis(hax(3,3),'image','off');
if ~any(isinterpolated(:)),
  res.hte_isinterpolated = text(nc/2,nr/2,'(None)','color','r',...
    'Parent',hax(3,3),'HorizontalAlignment','center');
end
imax(end+1) = hax(3,3);
bkgd_diagnostics.isinterpolated = isinterpolated;

%% fraction of pixels that were interpolated
bkgd_diagnostics.fracpx_interpolated = nnz(isinterpolated) / (nc*nr);

%% histogram of background center intensities for interpolated pixels

if any(isinterpolated(:)),
  counts = histc(ann.background_center(isinterpolated)',edges_bkgdcenter);
  counts = [counts(1:end-2),counts(end-1)+counts(end)];
  frac = counts / nnz(isinterpolated);
else
  frac = zeros(1,params.nbins_bkgdcenter);
end
res.hline_histisinterp = ...
  semilogy_with_zeros(hax(4,3),centers_bkgdcenter,frac,...
  [edges_bkgdcenter(1),edges_bkgdcenter(end)],1,...
  'k.-',params.linestyleparams_histisinterp{:});
res.hti_histisinterp = title(hax(4,3),'Bkgd Interpolated Px Intensities');

bkgd_diagnostics.frac_bkgdinterpolated = frac;

%% Image of background center log-likelihood ratio of foreground to background
expbgfgmodelfile = fullfile(settingsdir,analysis_protocol,dataloc_params.expbgfgmodelmatfilestr);
if isempty(expbgfgmodelfile) || ~exist(expbgfgmodelfile,'file'),
  llr = zeros(nr,nc);
else
  model = load(expbgfgmodelfile);
  llr = ExpBGFGModel_lik(model,ann.background_center);  
end
res.him_llr = imagesc(llr,'Parent',hax(3,2),params.clim_llr);
res.hti_llr = title(hax(3,2),'LLR Fg over Bg');
axis(hax(3,2),'image','off');
res.hcb_llr = colorbar('peer',hax(3,2),'Location','East');
cbpos = get(res.hcb_llr,'Position');
axpos = get(hax(3,2),'Position');
cbpos(1) = axpos(1)+axpos(3)+cbpos(3);
set(res.hcb_llr,'Position',cbpos);
imax(end+1) = hax(3,2);

%% histogram of background center llr
[edges_llr,centers_llr] = get_bins(params.nbins_llr,params.lim_llr,true,true);
counts = histc(llr(:)',edges_llr);
counts = [counts(1:end-2),counts(end-1)+counts(end)];
frac = counts / numel(dev);
res.hline_histllr = ...
  semilogy_with_zeros(hax(4,2),centers_llr,frac,...
  [edges_llr(1),edges_llr(end)],1,...
  'k.-',params.linestyleparams_llr{:});
res.hti_histllr = title(hax(4,2),'LLR Fg over Bg Hist');

bkgd_diagnostics.frac_bkgdcenter_llr = frac;
bkgd_diagnostics.edges_bkgdcenter_llr = edges_llr;

%% mode, average, prctiles of background center llr

% bin with highest count
[~,bin] = max(frac);
bkgd_diagnostics.histmode_bkgdcenter_llr = centers_bkgdcenter(bin);

% mean, std, prctiles
bkgd_diagnostics.prctiles_bkgdcenter_llr = prctile(llr(:),params.prctiles_compute);
bkgd_diagnostics.mean_bkgdcenter_llr = nanmean(llr(:));
bkgd_diagnostics.std_bkgdcenter_llr = nanstd(llr(:),1);

drawnow;

%% the following diagnostics involve sampling some frames from the video

moviefile = fullfile(expdir,dataloc_params.moviefilestr);
[readframe,nframes,fid] = get_readframe_fcn(moviefile);
trxfile = fullfile(expdir,dataloc_params.trxfilestr);
load(trxfile,'trx');

[edges_llrfore,centers_llrfore] = get_bins(params.nbins_llrfore,params.lim_llrfore,true,true);
[edges_diffim,centers_diffim] = get_bins(params.nbins_diffim,params.lim_diffim,true,true);
[edges_imfore,centers_imfore] = get_bins(params.nbins_imfore,params.lim_imfore,true,true);

imfore_counts = zeros(1,params.nbins_imfore);
llrfore_counts = zeros(1,params.nbins_llrfore);
diffim_counts = zeros(1,params.nbins_diffim);

sampleframes = round(linspace(1,nframes,params.nframessample));
frac_fg = 0;
nfliessample = 0;
minllrperfly = [];
maxllrperfly = [];
meanllrperfly = [];
mean_imfore = 0;
std_imfore = 0;
mean_diffim = 0;
std_diffim = 0;
mean_llrfore = 0;
std_llrfore = 0;
for i = 1:params.nframessample,
  
  t = sampleframes(i);
  fprintf('Sampling frame %d (%d/%d)\n',t,i,params.nframessample);
  
  % current image
  im = double(readframe(t));

  % background subtraction
  [isfore,diffim] = BackSub(im,ann);
  % use GMM with the tracked fly positions as initialization
  [cc,nfliescurr] = AssignPixels(isfore,diffim,trx,t);
  nfliessample = nfliessample + nfliescurr;
  
  % experiment-wide llr
  if ~exist('model','var'),
    llr = zeros(nr,nc);
  else
    llr = ExpBGFGModel_lik(model,im);
  end
  
  %% mean, std of imfore
  mean_imfore = mean_imfore + mean(im(isfore));
  std_imfore = std_imfore + mean(im(isfore).^2);
  
  %% histogram of imfore
  countscurr = histc(im(isfore)',edges_imfore);
  countscurr = [countscurr(1:end-2),countscurr(end-1)+countscurr(end)];
  imfore_counts = imfore_counts + countscurr;
  
  %% mean, std of diffim
  mean_diffim = mean_diffim + mean(diffim(:));
  std_diffim = std_diffim + mean(diffim(:).^2);
  
  %% histogram of difference image
  countscurr = histc(diffim(:)',edges_diffim);
  countscurr = [countscurr(1:end-2),countscurr(end-1)+countscurr(end)];
  diffim_counts = diffim_counts + countscurr;
  
  %% number of foreground pixels
  nforecurr = nnz(isfore);
  frac_fg = frac_fg + nforecurr;
  
  %% histogram of log-likelihood ratio of foreground to background for some flies
  countscurr = histc(llr(isfore)',edges_llrfore);
  countscurr = [countscurr(1:end-2),countscurr(end-1)+countscurr(end)];
  llrfore_counts = llrfore_counts + countscurr;
  
  %% mean, std of llrfore
  mean_llrfore = mean_llrfore + mean(llr(isfore));
  std_llrfore = std_llrfore + mean(llr(isfore).^2);

  %% Average, min, max log-likelihood ratio for some flies
  for j = 1:nfliescurr,
    llrcurr = llr(cc==j)';
    minllrperfly(end+1) = min(llrcurr); %#ok<AGROW>
    maxllrperfly(end+1) = max(llrcurr); %#ok<AGROW>
    meanllrperfly(end+1) = mean(llrcurr); %#ok<AGROW>
  end
  
end

%% imfore summary
bkgd_diagnostics.frac_imfore = imfore_counts / params.nframessample / (nr*nc);
bkgd_diagnostics.edges_imfore = edges_imfore;
% mode
[~,bin] = max(imfore_counts);
bkgd_diagnostics.histmode_imfore = centers_imfore(bin);
% compute mean, std
mean_imfore = mean_imfore / params.nframessample;
std_imfore = sqrt(std_imfore / params.nframessample - mean_imfore^2);
bkgd_diagnostics.mean_imfore = mean_imfore;
bkgd_diagnostics.std_imfore = std_imfore;
% approximate percentiles
bkgd_diagnostics.prctiles_imfore = nan(1,numel(params.prctiles_compute));
z = cumsum(bkgd_diagnostics.frac_imfore)*100;
for i = 1:numel(params.prctiles_compute),
  bin = find(z>=params.prctiles_compute(i),1);
  if isempty(bin),
    bkgd_diagnostics.prctiles_imfore(i) = edges_imfore(end-1);
  else
    bkgd_diagnostics.prctiles_imfore(i) = centers_imfore(end);
  end
end

% plot histogram
res.hline_imfore = ...
  semilogy_with_zeros(hax(1,4),centers_imfore,bkgd_diagnostics.frac_imfore,...
  [centers_imfore(1),centers_imfore(end)],1,...
  'k.-',params.linestyleparams_histimfore{:});
res.hti_histimfore = title(hax(1,4),'Foreground Px Intensities');

%% diffim summary
bkgd_diagnostics.frac_diffim = diffim_counts / params.nframessample / (nr*nc);
bkgd_diagnostics.edges_diffim = edges_diffim;
% mode
[~,bin] = max(diffim_counts);
bkgd_diagnostics.histmode_diffim = centers_diffim(bin);
% compute mean, std
mean_diffim = mean_diffim / params.nframessample;
std_diffim = sqrt(std_diffim / params.nframessample - mean_diffim^2);
bkgd_diagnostics.mean_diffim = mean_diffim;
bkgd_diagnostics.std_diffim = std_diffim;
% approximate percentiles
bkgd_diagnostics.prctiles_diffim = nan(1,numel(params.prctiles_compute));
z = cumsum(bkgd_diagnostics.frac_diffim)*100;
for i = 1:numel(params.prctiles_compute),
  bin = find(z>=params.prctiles_compute(i),1);
  if isempty(bin),
    bkgd_diagnostics.prctiles_diffim(i) = edges_diffim(end-1);
  else
    bkgd_diagnostics.prctiles_diffim(i) = centers_diffim(end);
  end
end

% plot histogram
res.hline_diffim = ...
  semilogy_with_zeros(hax(2,4),centers_diffim,bkgd_diagnostics.frac_diffim,...
  [centers_diffim(1),centers_diffim(end)],1,...
  'k.-',params.linestyleparams_histdiffim{:});
res.hti_histdiffim = title(hax(2,4),'Difference Image');

%% llrfore summary
bkgd_diagnostics.frac_llrfore = llrfore_counts / frac_fg;
bkgd_diagnostics.edges_llrfore = edges_llrfore;
% mode
[~,bin] = max(llrfore_counts);
bkgd_diagnostics.histmode_llrfore = centers_llrfore(bin);
% compute mean, std
mean_llrfore = mean_llrfore / params.nframessample;
std_llrfore = sqrt(std_llrfore / params.nframessample - mean_llrfore^2);
bkgd_diagnostics.mean_llrfore = mean_llrfore;
bkgd_diagnostics.std_llrfore = std_llrfore;
% approximate percentiles
bkgd_diagnostics.prctiles_llrfore = nan(1,numel(params.prctiles_compute));
z = cumsum(bkgd_diagnostics.frac_llrfore)*100;
for i = 1:numel(params.prctiles_compute),
  bin = find(z>=params.prctiles_compute(i),1);
  if isempty(bin),
    bkgd_diagnostics.prctiles_llrfore(i) = edges_llrfore(end-1);
  else
    bkgd_diagnostics.prctiles_llrfore(i) = centers_llrfore(end);
  end
end

bkgd_diagnostics.frac_fg = frac_fg / params.nframessample / (nr*nc);
bkgd_diagnostics.frac_minllrperfly = hist(minllrperfly,centers_llrfore)/nfliessample;
bkgd_diagnostics.frac_maxllrperfly = hist(maxllrperfly,centers_llrfore)/nfliessample;
bkgd_diagnostics.frac_meanllrperfly = hist(meanllrperfly,centers_llrfore)/nfliessample;

% mode
[~,bin] = max(bkgd_diagnostics.frac_minllrperfly);
bkgd_diagnostics.histmode_minllrperfly = centers_llrfore(bin);
[~,bin] = max(bkgd_diagnostics.frac_maxllrperfly);
bkgd_diagnostics.histmode_maxllrperfly = centers_llrfore(bin);
[~,bin] = max(bkgd_diagnostics.frac_meanllrperfly);
bkgd_diagnostics.histmode_meanllrperfly = centers_llrfore(bin);

% compute mean
bkgd_diagnostics.mean_minllrperfly = mean(minllrperfly);
bkgd_diagnostics.mean_meanllrperfly = mean(meanllrperfly);
bkgd_diagnostics.mean_maxllrperfly = mean(maxllrperfly);

% plot histograms
llr_frac = [bkgd_diagnostics.frac_llrfore
  bkgd_diagnostics.frac_minllrperfly
  bkgd_diagnostics.frac_meanllrperfly
  bkgd_diagnostics.frac_maxllrperfly];
res.hline_llrfore = ...
  semilogy_with_zeros(hax(3,4),centers_llrfore,llr_frac,[edges_llrfore(1),edges_llrfore(end)],1,'.-');
% hold(hax(3,4),'on');
% res.hline_minllrperfly = ...
%   semilogy_with_zeros(hax(3,4),centers_llrfore,bkgd_diagnostics.frac_minllrperfly,...
%   [centers_llrfore(1),centers_llrfore(end)],1,...
%   '.-','color',colors(2,:));
% res.hline_meanllrperfly = ...
%   semilogy_with_zeros(hax(3,4),centers_llrfore,bkgd_diagnostics.frac_meanllrperfly,...
%   [centers_llrfore(1),centers_llrfore(end)],1,...
%   '.-','color',colors(3,:));
% res.hline_maxllrperfly = ...
%   semilogy_with_zeros(hax(3,4),centers_llrfore,bkgd_diagnostics.frac_maxllrperfly,...
%   [centers_llrfore(1),centers_llrfore(end)],1,...
%   '.-','color',colors(4,:));
res.hti_histllrfore = title(hax(3,4),'Foreground LLR');
res.hlegend_histllrfore = legend([res.hline_llrfore],{'All fg','Min per fly','Mean per fly','Max per fly'});

outer_axpos = get(hax(3,4),'OuterPosition');
axpos = get(hax(3,4),'Position');
legpos = get(res.hlegend_histllrfore,'Position');
legpos(1) = axpos(1)+axpos(3)-legpos(3);
legpos(2) = outer_axpos(2)-legpos(4);
set(res.hlegend_histllrfore,'Position',legpos);

%% link axes
delete(hax(4,4));
linkaxes(imax);
fclose(fid);

%% save image

savename = fullfile(expdir,dataloc_params.bkgddiagnosticsimagefilestr);
if exist(savename,'file'),
  delete(savename);
end
save2png(savename,hfig);
  
%% save to mat file
bkgddiagnosticsmatfilename = fullfile(expdir,dataloc_params.bkgddiagnosticsmatfilestr);
save(bkgddiagnosticsmatfilename,'-struct','bkgd_diagnostics');

%% write to text file

bkgddiagnosticsfilename = fullfile(expdir,dataloc_params.bkgddiagnosticsfilestr);
fid1 = fopen(bkgddiagnosticsfilename,'w');
fns = {'bkgdcenter','bkgddev','alwaysbkgd','bkgdcenter_llr',...
  'imfore','diffim','llrfore','minllrperfly','meanllrperfly','maxllrperfly'};
statfns = {'histmode','mean','std','prctiles','frac'};
for i = 1:numel(fns),
  for j = 1:numel(statfns),
    fn = fns{i};
    statfn = statfns{j};
    fn1 = sprintf('%s_%s',statfn,fn);
    if ~isfield(bkgd_diagnostics,fn1),
      %fprintf('No field %s\n',fn1);
      continue;
    end
    fprintf(fid1,'%s',fn1);
    fprintf(fid1,',%f',bkgd_diagnostics.(fn1));
    fprintf(fid1,'\n');
  end
end
fclose(fid1);

function [edges,centers] = get_bins(nbins,lim,catchless,catchmore)

nextra = double(catchless) + double(catchmore);
edges = linspace(lim(1),lim(2),nbins+1-nextra);
if catchless,
  edges = [-inf,edges];
end
if catchmore,
  edges = [edges,inf];
end
centers = (edges(1:end-1)+edges(2:end))/2;
if catchmore,
  centers(end) = 2*centers(end-1) - centers(end-2);
end
if catchless,
  centers(1) = 2*centers(2) - centers(3);
end
