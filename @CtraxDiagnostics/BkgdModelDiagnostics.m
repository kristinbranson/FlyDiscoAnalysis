function res = BkgdModelDiagnostics(obj,varargin)

bkgd_diagnostics = struct;

[hax,hfig,figpos,...
  expdir,n,...
  edges_bkgdcenter,nbins_bkgdcenter,lim_bkgdcenter,...
  linestyleparams_bkgdcenter,linestyleparams_histalwaysbkgd,...
  edges_bkgddev,nbins_bkgddev,lim_bkgddev,clim_bkgddev,...
  linestyleparams_bkgddev,...
  edges_fracframesisback,nbins_fracframesisback,lim_fracframesisback,...
  linestyleparams_histfracframesisback,...
  edges_llr,nbins_llr,lim_llr,clim_llr,...
  linestyleparams_llr,...
  linestyleparams_histisinterp] = ...
  myparse(varargin,...
  'hax',[],'hfig',[],'figpos',[],...
  'expdir',[],'n',[],...
  'edges_bkgdcenter',[],...
  'nbins_bkgdcenter',50,...
  'lim_bkgdcenter',[0,255],...
  'linestyleparams_bkgdcenter',{},...
  'linestyleparams_histalwaysbkgd',{},...
  'edges_bkgddev',[],...
  'nbins_bkgddev',21,...
  'lim_bkgddev',[0,20],...
  'clim_bkgddev',[0,20],...
  'linestyleparams_bkgddev',{},...
  'edges_fracframesisback',[],...
  'nbins_fracframesisback',50,...
  'lim_fracframesisback',[0,1],...
  'linestyleparams_histfracframesisback',{},...
  'edges_llr',[],...
  'nbins_llr',50,...
  'lim_llr',[-60,60],...
  'clim_llr',[-60,60],...
  'linestyleparams_llr',{},...
  'linestyleparams_histisinterp',{});

% choose experiments to plot
if isempty(n),
  if isempty(expdir),
    n = 1;
  else
    n = obj.expdir2n(expdir);
  end
end

nax_r = 2;
nax_c = 6;
axparams = {nax_r,nax_c,.05};
[hax1,hfig] = get_axes(hax,hfig,'axparams',axparams,'figpos',figpos);
res.hax = hax1; res.hfig = hfig;
hax = reshape(hax1(1:nax_r*nax_c),[nax_r,nax_c]);

imax = [];

%% image of background center
bkgd_diagnostics.background_center = obj.anns{n}.background_center;
res.him_bkgdcenter = image(repmat(uint8(obj.anns{n}.background_center),[1,1,3]),'Parent',hax(1,1));
res.hti_bkgdcenter = title(hax(1,1),'Background Center');
axis(hax(1,1),'image','xy','off');
imax(end+1) = hax(1,1);

%% histogram of background center intensity
if isempty(edges_bkgdcenter),
  edges_bkgdcenter = linspace(lim_bkgdcenter(1),lim_bkgdcenter(2),nbins_bkgdcenter+1);
else
  nbins_bkgdcenter = length(edges_bkgdcenter) - 1;
end
centers_bkgdcenter = (edges_bkgdcenter(1:end-1)+edges_bkgdcenter(2:end))/2;
counts = histc(obj.anns{n}.background_center(:)',edges_bkgdcenter);
counts = [counts(1:end-2),counts(end-1)+counts(end)];
frac = counts / numel(obj.anns{n}.background_center);

res.hline_histbkgdcenter = ...
  semilogy_with_zeros(hax(2,1),centers_bkgdcenter,frac,...
  [edges_bkgdcenter(1),edges_bkgdcenter(end)],1,...
  'k.-',linestyleparams_bkgdcenter{:});
res.hti_histbkgdcenter = title(hax(2,1),'Bkgd Ctr Px Intensities');
bkgd_diagnostics.frac_bkgdcenter = frac;
bkgd_diagnostics.edges_bkgdcenter = edges_bkgdcenter;

%% image of background deviation
dev = [];
devtype = '';
if isfield(obj.anns{n},'background_mad'),
  dev = obj.anns{n}.background_mad;
  devtype = 'MAD';
elseif isfield(obj.anns{n},'background_std'),
  dev = obj.anns{n}.background_std;
  devtype = 'Std';
elseif isfield(obj.anns{n},'background_dev'),
  dev = obj.anns{n}.background_dev;
  devtype = 'Dev';
end

if ~isempty(dev),
  if isnan(clim_bkgddev(2)),
    clim_bkgddev(2) = max(dev(:));
  end
  if isnan(clim_bkgddev(1)),
    clim_bkgddev(1) = 0;
  end
  res.him_bkgddev = imagesc(dev,'Parent',hax(1,2),clim_bkgddev);
  res.hti_bkgddev = title(hax(1,2),['Background ',devtype]);
  res.hcb_bkgddev = colorbar('peer',hax(1,2),'Location','East','XColor','w','YColor','w');
  axis(hax(1,2),'image','xy','off');
  imax(end+1) = hax(1,2);
  
  %% histogram of background deviation
  if isempty(edges_bkgddev),
    edges_bkgddev = [linspace(lim_bkgddev(1),lim_bkgddev(2),nbins_bkgddev),inf];
  else
    %nbins_bkgddev = length(edges_bkgddev) - 1;
  end
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
  'k.-',linestyleparams_bkgddev{:});
  res.hti_histbkgddev = title(hax(2,2),['Bkgd ',devtype,' Hist']);

  bkgd_diagnostics.frac_bkgddev = frac;
  bkgd_diagnostics.edges_bkgddev = edges_bkgddev;

  %% maximum, percentiles of background median absolute deviation
  bkgd_diagnostics.max_bkgddev = max(dev(:));
  if isempty(prctiles_bkgddev_compute),
    bkgd_diagnostics.prctiles_bkgddev = [];
  else
    bkgd_diagnostics.prctiles_bkgddev = prctile(dev(:),prctiles_bkgddev_compute);
  end
  
else
  
  bkgd_diagnostics.frac_bkgddev = [];
  bkgd_diagnostics.edges_bkgddev = [];
  bkgd_diagnostics.max_bkgddev = nan;
  bkgd_diagnostics.prctiles_bkgddev = nan(size(prctiles_bkgddev_compute));
  
end

%% Image of background center overlaid with areas that are not always
% background (according to the alwaysbkgd mask) in red.
if isfield(obj.anns{n},'isarena'),
  isarena = obj.anns{n}.isarena~=0;
else
  isarena = true([obj.nrs(n),obj.ncs(n)]);
  if all(isfield(obj.anns{n},{'max_nonarena','background_center'})),
    isarena = isarena & (obj.anns{n}.background_center > obj.anns{n}.max_nonarena);
  end
  if all(isfield(obj.anns{n},{'min_nonarena','background_center'})),
    isarena = isarena & (obj.anns{n}.background_center < obj.anns{n}.min_nonarena);
  end
  if all(isfield(obj.anns{n},{'do_set_circular_arena',...
      'arena_center_x','arena_center_y','arena_radius'})) && ...
      obj.anns{n}.do_set_circular_arena ~= 0,
    [gridx,gridy] = meshgrid(0:obj.ncs(n)-1,0:obj.nrs(n)-1);
    isarena = isarena & ...
      ( (gridx - obj.anns{n}.arena_center_x).^2 + ...
      (gridy - obj.anns{n}.arena_center_y).^2 ) <= obj.anns{n}.arena_radius^2;
  end
end
tmp = obj.anns{n}.background_center;
tmp(~isarena) = tmp(~isarena)*2 + .5;
imalwaysbkgd = uint8(cat(3,tmp,repmat(obj.anns{n}.background_center,[1,1,2])));
res.him_alwaysbkgd = image(imalwaysbkgd,'Parent',hax(1,3));
res.hti_alwaysbkgd = title(hax(1,3),'Always Bkgd');
axis(hax(1,3),'image','xy','off');
if all(isarena),
  res.hte_alwaysbkgd = text(obj.ncs(n)/2,obj.nrs(n)/2,'(None)','color','r',...
    'Parent',hax(1,3),'HorizontalAlignment','center');
end
imax(end+1) = hax(1,3);

bkgd_diagnostics.isarena = isarena;
bkgd_diagnostics.imalwaysbkgd = imalwaysbkgd;

%% fraction of pixels that are always background
bkgd_diagnostics.fracpx_alwaysbkgd = nnz(isarena) / (obj.nrs(n)*obj.ncs(n));

%% histogram of background center intensities for pixels that are always
% background
if any(~isarena),
  counts = histc(obj.anns{n}.background_center(~isarena)',edges_bkgdcenter);
  counts = [counts(1:end-2),counts(end-1)+counts(end)];
  frac = counts / nnz(~isarena);
else
  frac = zeros(1,nbins_bkgdcenter);
end
res.hline_histalwaysbkgd = ...
  semilogy_with_zeros(hax(2,3),centers_bkgdcenter,frac,...
  [edges_bkgdcenter(1),edges_bkgdcenter(end)],1,...
  'k.-',linestyleparams_histalwaysbkgd{:});
res.hti_histalwaysbkgd = title(hax(2,3),'Always Bkgd Px Intensities');

bkgd_diagnostics.frac_alwaysbkgd = frac;

%% image of fraction of frames that are used to estimate background model
if ~isfield(obj.anns{n},'use_expbgfgmodel') || ...
    obj.anns{n}.use_expbgfgmodel == 0 || ...
    ~isfield(obj.anns{n},'fracframesisback'),...
  fracframesisback = ones(obj.nrs(n),obj.ncs(n));
else
  fracframesisback = obj.anns{n}.fracframesisback;
end
res.him_fracframesisback = imagesc(fracframesisback,'Parent',hax(1,4),[0,1]);
res.hti_alwaysbkgd = title(hax(1,4),'Frac. Frames Used in Bkgd Model');
axis(hax(1,4),'image','xy','off');
res.hcb_fracframesisback = colorbar('peer',hax(1,4),'Location','East','XColor','w','YColor','w');
imax(end+1) = hax(1,4);

bkgd_diagnostics.fracframesisback = fracframesisback;

%% histogram of fraction of frames used to estimate background
if isempty(edges_fracframesisback),
  edges_fracframesisback = linspace(lim_fracframesisback(1),lim_fracframesisback(2),nbins_fracframesisback+1);
else
  nbins_fracframesisback = length(edges_fracframesisback) - 1;
end
centers_fracframesisback = (edges_fracframesisback(1:end-1)+edges_fracframesisback(2:end))/2;
counts = histc(fracframesisback(:)',edges_fracframesisback);
counts = [counts(1:end-2),counts(end-1)+counts(end)];
frac = counts / numel(fracframesisback);
res.hline_histfracframesisback = ...
  semilogy_with_zeros(hax(2,4),centers_fracframesisback,frac,...
  [edges_fracframesisback(1),edges_fracframesisback(end)],1,...
  'k.-',linestyleparams_histfracframesisback{:});
res.hti_histfracframesisback = title(hax(2,4),'Frac. Frames Used in Bkgd Model Hist');

bkgd_diagnostics.frac_fracframesisback = frac;
bkgd_diagnostics.edges_fracframesisback = edges_fracframesisback;

%% mean, min, max fraction of frames used to estimate background model
bkgd_diagnostics.mean_fracframesisback = mean(fracframesisback(:));
bkgd_diagnostics.min_fracframesisback = min(fracframesisback(:));
bkgd_diagnostics.max_fracframesisback = max(fracframesisback(:));

%% Image of background center overlaid with areas that were interpolated in red
if isfield(obj.anns{n},'min_frac_frames_isback'),
  isinterpolated = fracframesisback < obj.anns{n}.min_frac_frames_isback;
else
  isinterpolated = false(obj.nrs(n),obj.ncs(n));
end
b = bwboundaries(imdilate(isinterpolated,ones(2)),8);
res.him_isinterpolated = image(repmat(uint8(obj.anns{n}.background_center),[1,1,3]),'Parent',hax(1,6));
hold(hax(1,6),'on');
res.hli_isinterpolated = nan(size(b));
for i = 1:numel(b),
  res.hli_isinterpolated(i) = plot(hax(1,6),b{i}(:,1),b{i}(:,2),'r-');
end
res.hti_isinterpolated = title(hax(1,6),'Bkgd Interpolated');
axis(hax(1,6),'image','xy','off');
if ~any(isinterpolated),
  res.hte_isinterpolated = text(obj.ncs(n)/2,obj.nrs(n)/2,'(None)','color','r',...
    'Parent',hax(1,6),'HorizontalAlignment','center');
end
imax(end+1) = hax(1,6);
bkgd_diagnostics.isinterpolated = isinterpolated;

%% fraction of pixels that were interpolated
bkgd_diagnostics.fracpx_interpolated = nnz(isinterpolated) / (obj.ncs(n)*obj.nrs(n));

%% histogram of background center intensities for interpolated pixels

if any(isinterpolated),
  counts = histc(obj.anns{n}.background_center(isinterpolated)',edges_bkgdcenter);
  counts = [counts(1:end-2),counts(end-1)+counts(end)];
  frac = counts / nnz(isinterpolated);
else
  frac = zeros(1,nbins_bkgdcenter);
end
res.hline_histisinterp = ...
  semilogy_with_zeros(hax(2,6),centers_bkgdcenter,frac,...
  [edges_bkgdcenter(1),edges_bkgdcenter(end)],1,...
  'k.-',linestyleparams_histisinterp{:});
res.hti_histisinterp = title(hax(2,6),'Bkgd Interpolated Px Intensities');

linkaxes(imax);
bkgd_diagnostics.frac_bkgdinterpolated = frac;

%% Image of background center log-likelihood ratio of foreground to background
if isempty(obj.ExpBGFGModelMatFile) || ~exist(obj.ExpBGFGModelMatFile,'file'),
  llr = zeros(obj.nrs(n),obj.ncs(n));
else
  model = load(obj.ExpBGFGModelMatFile);
  llr = ExpBGFGModel_lik(model,obj.anns{n}.background_center);  
end
res.him_llr = imagesc(llr,'Parent',hax(1,5),clim_llr);
res.hti_llr = title(hax(1,5),'LLR Fg over Bg');
axis(hax(1,5),'image','xy','off');
res.hcb_llr = colorbar('peer',hax(1,5),'Location','East','XColor','w','YColor','w');
imax(end+1) = hax(1,5);

%% histogram of background center llr
if isempty(edges_llr),
  edges_llr = [-inf,linspace(lim_llr(1),lim_llr(2),nbins_llr-1),inf];
else
  %nbins_llr = length(edges_llr) - 1;
end
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
  'k.-',linestyleparams_llr{:});
res.hti_histllr = title(hax(2,5),'LLR Fg over Bg Hist');

bkgd_diagnostics.frac_bkgdcenter_llr = frac;
bkgd_diagnostics.edges_bkgdcenter_llr = edges_llr;

%% average, min, max of background center llr
bkgd_diagnostics.mean_bkgdcenter_llr = nanmean(llr(:));
bkgd_diagnostics.min_bkgdcenter_llr = min(llr(:));
bkgd_diagnostics.max_bkgdcenter_llr = max(llr(:));