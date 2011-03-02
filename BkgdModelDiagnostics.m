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
bkgdmodeldiagnosticsparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.bkgdmodeldiagnosticsparamsfilestr);
params = ReadParams(bkgdmodeldiagnosticsparamsfile);

%% set up figure

res = struct;
hfig = 1;
figure(hfig);
clf(hfig,'reset');
set(hfig,'Units','Pixels','Position',params.figpos);
nax_r = 2;
nax_c = 6;
hax = createsubplots(nax_r,nax_c,.05,hfig);
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
axis(hax(1,1),'image','xy','off');
imax(end+1) = hax(1,1);

%% histogram of background center intensity
edges_bkgdcenter = linspace(params.lim_bkgdcenter(1),params.lim_bkgdcenter(2),params.nbins_bkgdcenter+1);
centers_bkgdcenter = (edges_bkgdcenter(1:end-1)+edges_bkgdcenter(2:end))/2;
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
  axis(hax(1,2),'image','xy','off');
  imax(end+1) = hax(1,2);
  cbpos = get(res.hcb_bkgddev,'Position');
  axpos = get(hax(1,2),'Position');
  cbpos(1) = axpos(1)+axpos(3)+cbpos(3);
  set(res.hcb_bkgddev,'Position',cbpos);
  
  %% histogram of background deviation
  edges_bkgddev = [linspace(params.lim_bkgddev(1),params.lim_bkgddev(2),params.nbins_bkgddev),inf];
  if isinf(edges_bkgddev(end)),
    centers_bkgddev = (edges_bkgddev(1:end-2)+edges_bkgddev(2:end-1))/2;
    centers_bkgddev(end+1) = 2*centers_bkgddev(end) - centers_bkgddev(end-1);
  else
    centers_bkgddev = (edges_bkgddev(1:end-1)+edges_bkgddev(2:end))/2;
  end
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

  %% maximum, percentiles of background median absolute deviation
  bkgd_diagnostics.max_bkgddev = max(dev(:));
  if isempty(params.prctiles_bkgddev_compute),
    bkgd_diagnostics.prctiles_bkgddev = [];
  else
    bkgd_diagnostics.prctiles_bkgddev = prctile(dev(:),params.prctiles_bkgddev_compute);
  end
  
else
  
  bkgd_diagnostics.frac_bkgddev = [];
  bkgd_diagnostics.edges_bkgddev = [];
  bkgd_diagnostics.max_bkgddev = nan;
  bkgd_diagnostics.prctiles_bkgddev = nan(size(prctiles_bkgddev_compute));
  
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
tmp(~isarena) = tmp(~isarena)*2 + .5;
imalwaysbkgd = uint8(cat(3,tmp,repmat(ann.background_center,[1,1,2])));
res.him_alwaysbkgd = image(imalwaysbkgd,'Parent',hax(1,3));
res.hti_alwaysbkgd = title(hax(1,3),'Always Bkgd');
axis(hax(1,3),'image','xy','off');
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

%% image of fraction of frames that are used to estimate background model
if ~isfield(ann,'use_expbgfgmodel') || ...
    ann.use_expbgfgmodel == 0 || ...
    ~isfield(ann,'fracframesisback'),
  fracframesisback = ones(nr,nc);
else
  fracframesisback = ann.fracframesisback;
end
res.him_fracframesisback = imagesc(fracframesisback,'Parent',hax(1,4),[0,1]);
res.hti_alwaysbkgd = title(hax(1,4),'Frac. Frames Used in Bkgd Model');
axis(hax(1,4),'image','xy','off');
res.hcb_fracframesisback = colorbar('peer',hax(1,4),'Location','East');
cbpos = get(res.hcb_fracframesisback,'Position');
axpos = get(hax(1,4),'Position');
cbpos(1) = axpos(1)+axpos(3)+cbpos(3);
set(res.hcb_fracframesisback,'Position',cbpos);
imax(end+1) = hax(1,4);

bkgd_diagnostics.fracframesisback = fracframesisback;

%% histogram of fraction of frames used to estimate background
edges_fracframesisback = linspace(params.lim_fracframesisback(1),...
  params.lim_fracframesisback(2),params.nbins_fracframesisback+1);
centers_fracframesisback = (edges_fracframesisback(1:end-1)+edges_fracframesisback(2:end))/2;
counts = histc(fracframesisback(:)',edges_fracframesisback);
counts = [counts(1:end-2),counts(end-1)+counts(end)];
frac = counts / numel(fracframesisback);
res.hline_histfracframesisback = ...
  semilogy_with_zeros(hax(2,4),centers_fracframesisback,frac,...
  [edges_fracframesisback(1),edges_fracframesisback(end)],1,...
  'k.-',params.linestyleparams_histfracframesisback{:});
res.hti_histfracframesisback = title(hax(2,4),'Frac. Frames Used in Bkgd Model Hist');

bkgd_diagnostics.frac_fracframesisback = frac;
bkgd_diagnostics.edges_fracframesisback = edges_fracframesisback;

%% mean, min, max fraction of frames used to estimate background model
bkgd_diagnostics.mean_fracframesisback = mean(fracframesisback(:));
bkgd_diagnostics.min_fracframesisback = min(fracframesisback(:));
bkgd_diagnostics.max_fracframesisback = max(fracframesisback(:));

%% Image of background center overlaid with areas that were interpolated in red
if isfield(ann,'min_frac_frames_isback'),
  isinterpolated = fracframesisback < ann.min_frac_frames_isback;
else
  isinterpolated = false(nr,nc);
end
b = bwboundaries(imdilate(isinterpolated,ones(2)),8);
res.him_isinterpolated = image(repmat(uint8(ann.background_center),[1,1,3]),'Parent',hax(1,6));
hold(hax(1,6),'on');
res.hli_isinterpolated = nan(size(b));
for i = 1:numel(b),
  res.hli_isinterpolated(i) = plot(hax(1,6),b{i}(:,1),b{i}(:,2),'r-');
end

res.hti_isinterpolated = title(hax(1,6),'Bkgd Interpolated');
axis(hax(1,6),'image','xy','off');
if ~any(isinterpolated(:)),
  res.hte_isinterpolated = text(nc/2,nr/2,'(None)','color','r',...
    'Parent',hax(1,6),'HorizontalAlignment','center');
end
imax(end+1) = hax(1,6);
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
  semilogy_with_zeros(hax(2,6),centers_bkgdcenter,frac,...
  [edges_bkgdcenter(1),edges_bkgdcenter(end)],1,...
  'k.-',params.linestyleparams_histisinterp{:});
res.hti_histisinterp = title(hax(2,6),'Bkgd Interpolated Px Intensities');

bkgd_diagnostics.frac_bkgdinterpolated = frac;

%% Image of background center log-likelihood ratio of foreground to background
expbgfgmodelfile = fullfile(settingsdir,analysis_protocol,dataloc_params.expbgfgmodelmatfilestr);
if isempty(expbgfgmodelfile) || ~exist(expbgfgmodelfile,'file'),
  llr = zeros(nr,nc);
else
  model = load(expbgfgmodelfile);
  llr = ExpBGFGModel_lik(model,ann.background_center);  
end
res.him_llr = imagesc(llr,'Parent',hax(1,5),params.clim_llr);
res.hti_llr = title(hax(1,5),'LLR Fg over Bg');
axis(hax(1,5),'image','xy','off');
res.hcb_llr = colorbar('peer',hax(1,5),'Location','East');
cbpos = get(res.hcb_llr,'Position');
axpos = get(hax(1,5),'Position');
cbpos(1) = axpos(1)+axpos(3)+cbpos(3);
set(res.hcb_llr,'Position',cbpos);
imax(end+1) = hax(1,5);

%% histogram of background center llr
edges_llr = [-inf,linspace(params.lim_llr(1),params.lim_llr(2),params.nbins_llr-1),inf];
centers_llr = (edges_llr(2:end-2)+edges_llr(3:end-1))/2;
if isinf(edges_llr(end)),
  centerend = 2*centers_llr(end) - centers_llr(end-1);
else
  centerend = (edges_llr(end-1)+edges_llr(end))/2;
end
if isinf(edges_llr(1)),
  centerstart = 2*centers_llr(1) - centers_llr(2);
else
  centerstart = (edges_llr(1)+edges_llr(2))/2;
end
centers_llr = [centerstart,centers_llr,centerend];
counts = histc(llr(:)',edges_llr);
counts = [counts(1:end-2),counts(end-1)+counts(end)];
frac = counts / numel(dev);
res.hline_histllr = ...
  semilogy_with_zeros(hax(2,5),centers_llr,frac,...
  [edges_llr(1),edges_llr(end)],1,...
  'k.-',params.linestyleparams_llr{:});
res.hti_histllr = title(hax(2,5),'LLR Fg over Bg Hist');

bkgd_diagnostics.frac_bkgdcenter_llr = frac;
bkgd_diagnostics.edges_bkgdcenter_llr = edges_llr;

%% average, min, max of background center llr
bkgd_diagnostics.mean_bkgdcenter_llr = nanmean(llr(:));
bkgd_diagnostics.min_bkgdcenter_llr = min(llr(:));
bkgd_diagnostics.max_bkgdcenter_llr = max(llr(:));

%% the following diagnostics involve sampling some frames from the video

moviefile = fullfile(expdir,dataloc_params.moviefilestr);
[readframe,nframes,fid] = get_readframe_fcn(moviefile);
trxfile = fullfile(expdir,dataloc_params.trxfilestr);
load(trxfile,'trx');

edges_llrfore = [-inf,linspace(params.lim_llrfore(1),params.lim_llrfore(2),...
  params.nbins_llrfore-1),inf];
centers_llrfore = (edges_llrfore(2:end-2)+edges_llrfore(3:end-1))/2;
if isinf(edges_llrfore(end)),
  centerend = 2*centers_llrfore(end) - centers_llrfore(end-1);
else
  centerend = (edges_llrfore(end-1)+edges_llrfore(end))/2;
end
if isinf(edges_llrfore(1)),
  centerstart = 2*centers_llrfore(1) - centers_llrfore(2);
else
  centerstart = (edges_llrfore(1)+edges_llrfore(2))/2;
end
centers_llrfore = [centerstart,centers_llrfore,centerend];
llrfore_counts = zeros(1,params.nbins_llrfore);

sampleframes = round(linspace(1,nframes,params.nframessample));
diffim_counts = zeros(1,256);
centers = 0:255;
frac_fg = 0;
nfliessample = 0;
minllrperfly = [];
maxllrperfly = [];
meanllrperfly = [];
for i = 1:params.nframessample,
  
  t = sampleframes(i);
  
  % current image
  im = double(readframe(t));
%   
%   % current positions
%   trxcurr = obj.trx.getmovieidx(n).getframe(t);
%   nfliescurr = numel(trxcurr);
%   nfliessample = nfliessample + nfliescurr;
%   
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
  
  %% histogram of difference image
  countscurr = hist(diffim(:),centers);
  diffim_counts = diffim_counts + countscurr;
  
  %% number of foreground pixels
  nforecurr = nnz(isfore);
  frac_fg = frac_fg + nforecurr;
  
  %% histogram of log-likelihood ratio of foreground to background for some flies
  llrfore_counts = llrfore_counts + hist(llr(isfore)',centers_llrfore);
  
  %% Average, min, max log-likelihood ratio for some flies
  for j = 1:nfliescurr,
    llrcurr = llr(cc==j)';
    minllrperfly(end+1) = min(llrcurr); %#ok<AGROW>
    maxllrperfly(end+1) = max(llrcurr); %#ok<AGROW>
    meanllrperfly(end+1) = mean(llrcurr); %#ok<AGROW>
  end
  
end

bkgd_diagnostics.diffim_frac = diffim_counts / params.nframessample / (nr*nc);
bkgd_diagnostics.centers_diffim = centers;
bkgd_diagnostics.llrfore_frac = llrfore_counts / frac_fg;
bkgd_diagnostics.edges_llrfore = edges_llrfore;
bkgd_diagnostics.frac_fg = frac_fg / params.nframessample / (nr*nc);
bkgd_diagnostics.frac_minllrperfly = hist(minllrperfly,centers_llrfore)/nfliessample;
bkgd_diagnostics.frac_maxllrperfly = hist(maxllrperfly,centers_llrfore)/nfliessample;
bkgd_diagnostics.frac_meanllrperfly = hist(meanllrperfly,centers_llrfore)/nfliessample;

%% link axes
linkaxes(imax);
fclose(fid);


