%% set up path

%rootdatadir = '/groups/branson/bransonlab/projects/olympiad/WingTrackedData';
addpath ../misc;
addpath ../filehandling;
if ispc,
  addpath ../../Ctrax/trunk/matlab/netlab;
else
  addpath ../../Ctrax/matlab/netlab;
end

%addpath /groups/branson/home/bransonk/olympiad_bransonlab/FlyBowlAnalysis;
rootdatadir = '/groups/branson/home/bransonk/tracking/data/heberlein';
experiment_name = '11292012MVMM1';

%% parameters

moviefile = fullfile(rootdatadir,experiment_name,'movie.avi');
annfile = '';
%annfile = fullfile(rootdatadir,experiment_name,'movie.ufmf.ann');

% number of frames to sample when estimating the background model
bg_nframes = 200;

% frames to track
firstframetrack = 1;
lastframetrack = inf;

% background type
DARKBKGD = 1;
LIGHTBKGD = -1;
OTHERBKGD = 0;
bgmode = LIGHTBKGD;

% background subtraction threshold
bgthresh = .5;

% background subtraction normalization image
minbgnorm = 80;

% minimum area of a connected component in pixels
minccarea = 10;

% gmm em
gmmem_params.nrestarts_firstframe = 10;
gmmem_params.precision = .01;
gmmem_params.maxiters = 100;
gmmem_params.min_obsprior = .2;

% matching
err_params.dampen_priors = .1;
err_params.dampen_pos = .1;
err_params.dampen_theta = .8;
err_params.weightpos = 1;
err_params.weighttheta = 1;
err_params.weightarea = 1;

% resolving head-tail ambiguity
choose_orientations_params.velocity_angle_weight = .03;
choose_orientations_params.max_velocity_angle_weight = .18;

flycolors = {'r','b'};

DEBUG = 2;

%% read stuff in

[readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefile);
imheight = headerinfo.nr;
imwidth = headerinfo.nc;

%% compute background model

if ~isempty(annfile) && exist(annfile,'file'),
  bgmed = read_ann(annfile,'background_center');
  bgmed = uint8(round(reshape(bgmed,[imwidth,imheight])'));
else
  buffer = readframe(1);
  iscolor = false;
  if size(buffer,3) > 1,
    iscolor = true;
    buffer = rgb2gray(buffer);
  end
  buffer = repmat(buffer,[1,1,bg_nframes]);
  frames = round(linspace(1,nframes,bg_nframes));
  for i = 1:bg_nframes,
    t = frames(i);
    im = readframe(t);
    if iscolor,
      im = rgb2gray(im);
    end
    buffer(:,:,i) = im;
  end
  
  bgmed = uint8(median(single(buffer),3));
%   
%   
%   for i = 1:bg_nframes,
%     t = frames(i);
%     im = readframe(t);
%     if iscolor,
%       im = rgb2gray(im);
%     end
%     dbkgd = imabsdiff(im,bgmed);
%     buffer(:,:,i) = dbkgd;
%   end
%   
%   bgmad = uint8(median(single(buffer),3));
  
end


%% get regions of interest

hfig = 1;
figure(hfig);
clf;
imagesc(bgmed,[0,255]);
axis image;
colormap jet;
hax = gca;
title('Click to add next region of interest. Double-click on the polygon to continue.');

idxroi = zeros(size(bgmed));
idxtethered = zeros(size(bgmed));
rois = {};
tethered = {};
inrois = {};
intethered = {};
roibbs = zeros(0,4);
tetheredbbs = zeros(0,4);
[XGRID,YGRID] = meshgrid(1:imwidth,1:imheight);

bgmed_fixed = bgmed;

while true,
  
  title('Click to add next region of interest. Double-click on the polygon to continue.');
  h = impoly(hax,'Closed',true);
  position = wait(h);
  hold on;
  hplot = plot(position(:,1),position(:,2),'k.-');
  if ishandle(h),
    delete(h);
  end
  
  title('Click to outline tethered fly. Double-click on the polygon to continue.');
  h = impoly(hax,'Closed',true);
  tethered_position = wait(h);
  hold on;
  hplot_tethered = plot(tethered_position(:,1),tethered_position(:,2),'w.-');
  if ishandle(h),
    delete(h);
  end
  
  res = questdlg('Choose next action:','ROI next action','Add another ROI','Finish','Undo','Add another ROI');
  switch res,
    case {'Finish','Add another ROI'},
      rois{end+1} = position; %#ok<SAGROW>
      bb = [floor(max(1,floor(min(position(:,1))))),ceil(min(imwidth,(max(position(:,1))))),...
        floor(max(1,floor(min(position(:,2))))),ceil(min(imheight,(max(position(:,2)))))];
      roibbs(end+1,:) = bb; %#ok<SAGROW>
      in = inpolygon(XGRID,YGRID,position(:,1),position(:,2));
      idxroi(in) = numel(rois);
      inrois{end+1} = in(bb(3):bb(4),bb(1):bb(2)); %#ok<SAGROW>

      tethered{end+1} = tethered_position; %#ok<SAGROW>
      tetheredbbs(end+1,:) = bb; %#ok<SAGROW>
      in = inpolygon(XGRID,YGRID,tethered_position(:,1),tethered_position(:,2));
      idxtethered(in) = numel(tethered);
      intethered{end+1} = in(bb(3):bb(4),bb(1):bb(2)); %#ok<SAGROW>
      
      if strcmp(res,'Finish'),
        break;
      end
    case 'Undo',
      if ishandle(hplot),
        delete(hplot);
      end
  end
end

nrois = numel(rois);

save(fullfile(expdir,'roiinfo.mat'),'idxroi','idxtethered',...
  'rois','tethered','inrois','intethered','roibbs','tetheredbbs','nrois');

%% normalization image

inroi = false(size(bgmed));
for i = 1:nrois,
  inroi(roibbs(i,3):roibbs(i,4),roibbs(i,1):roibbs(i,2)) = ...
    inroi(roibbs(i,3):roibbs(i,4),roibbs(i,1):roibbs(i,2)) | inrois{i};
end

bgnorm = double(bgmed);
bgnorm(~inroi) = inf;
for i = 1:numel(rois),
  tmp = bgnorm(roibbs(i,3):roibbs(i,4),roibbs(i,1):roibbs(i,2));
  tmp2 = roifill(tmp,intethered{i});
  bgnorm(roibbs(i,3):roibbs(i,4),roibbs(i,1):roibbs(i,2)) = tmp2;
end
bgnorm = max(minbgnorm,bgnorm);


%% look at the background subtraction results

colors = jet(1000);
colors = colors(randperm(1000),:);
hb = [];
isfirst = true;

for t = randperm(nframes),

im = readframe(t);
if iscolor,
  im = rgb2gray(im);
end

switch bgmode,
  
  case DARKBKGD,
    dbkgd = imsubtract(im,bgmed);
  case LIGHTBKGD,
    dbkgd = imsubtract(bgmed,im);
  case OTHERBKGD,
    dbkgd = imabsdiff(im,bgmed);
end
dbkgd = double(dbkgd)./bgnorm;

imthresh = dbkgd >= bgthresh;
imthresh = bwareaopen(imthresh,minccarea);

b = bwboundaries(imthresh);
if isfirst || ~exist('him','var') || ~ishandle(him),
  clf;
  him = imagesc(im,[0,255]); axis image;
  colormap gray;
  hax = gca;
  hold on;
else
  set(him,'CData',im);
end

set(hb(ishandle(hb)),'XData',nan,'YData',nan);

for i = 1:numel(b),
  if isfirst || ~exist('h','var') || numel(hb) < i || ~ishandle(hb(i)),
    hb(i) = plot(hax,b{i}(:,2),b{i}(:,1),'-','Color',colors(i,:));
  else
    set(hb(i),'XData',b{i}(:,2),'YData',b{i}(:,1));
  end
end

isfirst = false;
input(num2str(t));

end

%% initialize

nframes_track = min(lastframetrack,nframes)-firstframetrack+1;
trxx = nan(2,nrois,nframes_track);
trxy = nan(2,nrois,nframes_track);
trxa = nan(2,nrois,nframes_track);
trxb = nan(2,nrois,nframes_track);
trxtheta = nan(2,nrois,nframes_track);
trxarea = nan(2,nrois,nframes_track);
% mixpred = gmm(2,2,'full');
% mixpred = repmat(mixpred,[1,nrois]);
% predx = nan(2,nrois);
% predy = nan(2,nrois);
% predtheta = nan(2,nrois);
% predarea = nan(2,nrois);
% istouching = nan(nrois,nframes_track);

hell = nan(2,nrois);
htrx = nan(2,nrois);

% gmmfit_nbadpriors = zeros(1,nrois);

%% find the tethered flies

for i = 1:nrois,
  
  im = double(bgmed(roibbs(i,3):roibbs(i,4),roibbs(i,1):roibbs(i,2)));
  bg = bgnorm(roibbs(i,3):roibbs(i,4),roibbs(i,1):roibbs(i,2));
  switch bgmode,
    case DARKBKGD,
      dbkgd = imsubtract(im,bg);
    case LIGHTBKGD,
      dbkgd = imsubtract(bg,im);
    case OTHERBKGD,
      dbkgd = imabsdiff(im,bg);
  end
  dbkgd(~intethered{i}) = 0;
  dbkgd = double(dbkgd) ./ bg;
  isfore = dbkgd >= bgthresh;
  %isfore = bwareaopen(isfore,minccarea);
  cc = bwconncomp(isfore);
  if cc.NumObjects > 1,
    areas = cellfun(@numel,cc.PixelIdxList);
    [area,j] = max(areas);
  else
    j = 1;
    area = numel(cc.PixelIdxList{1});
  end
  [y,x] = ind2sub(cc.ImageSize,cc.PixelIdxList{j});
  [mu,S] = weighted_mean_cov([x(:),y(:)],dbkgd(cc.PixelIdxList{j}));
  mu(1) = mu(1);
  mu(2) = mu(2);
  trxx(2,i,:) = mu(1);
  trxy(2,i,:) = mu(2);
  [a,b,theta] = cov2ell(S);
  trxa(2,i,:) = a;
  trxb(2,i,:) = b;
  trxtheta(2,i,:) = theta;
  trxarea(2,i,:) = area;
  
end

%% resolve head/tail
clf;
for i = 1:nrois,

  hold off;
  im = bgmed(roibbs(i,3):roibbs(i,4),roibbs(i,1):roibbs(i,2));
  imagesc(im,[0,255]);
  colormap jet;
  axis image;
  hold on;
  drawflyo(trxx(2,i,1),trxy(2,i,1),trxtheta(2,i,1),trxa(2,i,1)/2,trxb(2,i,1)/2,'k');
  %drawflyo(trxx(2,i,1)-roibbs(i,1)+1,trxy(2,i,1)-roibbs(i,3)+1,trxtheta(2,i,1),trxa(2,i,1)/2,trxb(2,i,1)/2,'k');

    
  res = questdlg('Flip fly?');
  if strcmpi(res,'Yes'),
    trxtheta(2,i,:) = modrange(trxtheta(2,i,:) + pi,-pi,pi);
  elseif strcmpi(res,'Cancel'),
    break;
  end
  
end


%%

for t = t:min(lastframetrack,nframes),
  
  %% backsub
  
  iframe = t - firstframetrack+1;
  
  % read in frame
  if DEBUG >= 2,
    fprintf('Frame %d\n',t);
    tic;
  end
  im = readframe(t);
  if iscolor,
    im = rgb2gray(im);
  end
  if DEBUG >= 2,
    fprintf('Reading frame %d took %f seconds\n',t,toc);
    tic;
  end
  
  % subtract off background
  switch bgmode,
    case DARKBKGD,
      dbkgd = imsubtract(im,bgmed);
    case LIGHTBKGD,
      dbkgd = imsubtract(bgmed,im);
    case OTHERBKGD,
      dbkgd = imabsdiff(im,bgmed);
  end
  
  dbkgd = double(dbkgd) ./ bgnorm;
  
  % threshold
  isfore = dbkgd >= bgthresh;
  
  if DEBUG >= 2,
    fprintf('Backsub took %f seconds\n',toc);
  end

  
  if DEBUG,
    if t == firstframetrack || ~exist('him','var') || ~ishandle(him),
      hold off;
      him = imagesc(im,[0,255]);
      axis image;
      colormap gray;
      hold on;
      isnewplot = true;
    else
      set(him,'CData',im);
      isnewplot = false;
    end
    title(num2str(t));
  end
 
  %% loop over rois
  
  for roii = 1:nrois,
    
    if DEBUG >= 2,
      fprintf('Frame %d, roii %d\n',t,roii);
    end
    
    %% get bounding box & roi
    
    if DEBUG >= 2,
      tic;
    end
    
    roibb = roibbs(roii,:);
    offx = roibb(1)-1;
    offy = roibb(3)-1;
    isforebb = isfore(roibb(3):roibb(4),roibb(1):roibb(2));
    dbkgdbb = double(dbkgd(roibb(3):roibb(4),roibb(1):roibb(2)));
    dbkgdbb(~inrois{roii}) = 0;
    isforebb(~inrois{roii}) = false;
    
    %% fit connected components
    cc = bwconncomp(isforebb);
    
    cc.Area = cellfun(@numel,cc.PixelIdxList);
    [~,j] = max(cc.Area);
    ncc_bigenough = nnz(cc.Area >= minccarea);

    if DEBUG >= 2,
      fprintf('Connected components analysis took %f seconds\n',toc);
      tic;
    end

    
    %% fit gaussian

    if cc.NumObjects == 0,
      trxx(1,roii,iframe) = nan;
      trxy(1,roii,iframe) = nan;
      trxa(1,roii,iframe) = 0;
      trxb(1,roii,iframe) = 0;
      trxtheta(1,roii,iframe) = nan;
      trxarea(1,roii,iframe) = 0;
      fprintf('No fly detected in frame %d, roi %d\n',t,roii);
    else
      [y,x] = ind2sub(size(isforebb),cc.PixelIdxList{j});
      w = dbkgdbb(cc.PixelIdxList{j});
      [mu,S] = weighted_mean_cov([x(:),y(:)],w);
      [a,b,theta] = cov2ell(S);
      trxx(1,roii,iframe) = mu(1);
      trxy(1,roii,iframe) = mu(2);
      trxa(1,roii,iframe) = a;
      trxb(1,roii,iframe) = b;
      trxtheta(1,roii,iframe) = theta;
      trxarea(1,roii,iframe) = cc.Area(j);
    end

    %% plot
    
    if DEBUG,
      for i = 1:2,
        if isnewplot || ~exist('hell','var') || ~ishandle(hell(i,roii)),
          hell(i,roii) = drawellipse(trxx(i,roii,iframe)+offx,trxy(i,roii,iframe)+offy,trxtheta(i,roii,iframe),trxa(i,roii,iframe),trxb(i,roii,iframe),[flycolors{i},'-']);
        else
          updateellipse(hell(i,roii),trxx(i,roii,iframe)+offx,trxy(i,roii,iframe)+offy,trxtheta(i,roii,iframe),trxa(i,roii,iframe),trxb(i,roii,iframe));
        end
        if isnewplot || ~exist('htrx','var') || ~ishandle(htrx(i,roii)),
          htrx(i,roii) = plot(squeeze(trxx(i,roii,max(iframe-30,1):iframe)+offx),squeeze(trxy(i,roii,max(iframe-30,1):iframe)+offy),[flycolors{i},'.-']);
        else
          set(htrx(i,roii),'XData',squeeze(trxx(i,roii,max(iframe-30,1):iframe)+offx),...
            'YData',squeeze(trxy(i,roii,max(iframe-30,1):iframe)+offy));
        end
      end
      %input('');
    end
    
  end
  
  fprintf('Done with frame %d\n',t);
  if DEBUG,
    drawnow;
  end      
  
end

%% correct for bounding box of rois
for roii = 1:nrois,
  roibb = roibbs(roii,:);
  trxx(:,roii,:) = trxx(:,roii,:) + roibb(1) - 1;
  trxy(:,roii,:) = trxy(:,roii,:) + roibb(3) - 1;
end

%% resolve head/tail ambiguity

for roii = 1:nrois,

  i = 1;
  x = reshape(trxx(i,roii,:),[1,nframes_track]);
  y = reshape(trxy(i,roii,:),[1,nframes_track]);
  theta = reshape(trxtheta(i,roii,:),[1,nframes_track]);
  trxtheta(i,roii,:) = choose_orientations(x,y,theta,choose_orientations_params.velocity_angle_weight,choose_orientations_params.max_velocity_angle_weight);
  
end

%% plot updated results

for iframe = 1:nframes_track,
  
  t = iframe+firstframetrack-1;
  im = readframe(t);

  if t == firstframetrack || ~exist('him','var') || ~ishandle(him),
    hold off;
    him = imagesc(im,[0,255]);
    axis image;
    colormap gray;
    hold on;
    isnewplot = true;
  else
    set(him,'CData',im);
    isnewplot = false;
  end
  title(num2str(t));
  
  for roii = 1:nrois,
    for i = 1:2,
      if isnewplot || ~exist('hell','var') || ~ishandle(hell(i,roii)),
        hell(i,roii) = drawflyo(trxx(i,roii,iframe),trxy(i,roii,iframe),trxtheta(i,roii,iframe),trxa(i,roii,iframe)/2,trxb(i,roii,iframe)/2,[flycolors{i},'-']);
      else
        updatefly(hell(i,roii),trxx(i,roii,iframe),trxy(i,roii,iframe),trxtheta(i,roii,iframe),trxa(i,roii,iframe)/2,trxb(i,roii,iframe)/2);
      end
      if isnewplot || ~exist('htrx','var') || ~ishandle(htrx(i,roii)),
        htrx(i,roii) = plot(squeeze(trxx(i,roii,max(iframe-30,1):iframe)),squeeze(trxy(i,roii,max(iframe-30,1):iframe)),[flycolors{i},'.-']);
      else
        set(htrx(i,roii),'XData',squeeze(trxx(i,roii,max(iframe-30,1):iframe)),...
          'YData',squeeze(trxy(i,roii,max(iframe-30,1):iframe)));
      end
    end
    %input('');
  end
  drawnow;
  
end

%% save results

[~,~,ext] = fileparts(moviefile);
if strcmpi(ext,'.ufmf'),
  timestamps = ufmf_read_timestamps(headerinfo,1,nframes);
else
  timestamps = (0:nframes-1)/ headerinfo.FrameRate;
end
timestamp_tracked = now;
 
% trx = struct('x',mat2cell(reshape(trxx,[2*nrois,nframes_track]),ones(1,2*nrois),nframes_track),...
%   'y',mat2cell(reshape(trxy,[2*nrois,nframes_track]),ones(1,2*nrois),nframes_track),...
%   'a',mat2cell(reshape(trxa/2,[2*nrois,nframes_track]),ones(1,2*nrois),nframes_track),...
%   'b',mat2cell(reshape(trxb/2,[2*nrois,nframes_track]),ones(1,2*nrois),nframes_track),...
%   'theta',mat2cell(reshape(trxtheta,[2*nrois,nframes_track]),ones(1,2*nrois),nframes_track));
trx = struct;
for i = 1:nrois,
  for jj = 1:2,
    j = sub2ind([2,nrois],jj,i);
    trx(j).x = reshape(trxx(jj,i,:),[1,nframes_track]);
    trx(j).y = reshape(trxy(jj,i,:),[1,nframes_track]);
    % saved quarter-major, quarter-minor axis
    trx(j).a = reshape(trxa(jj,i,:)/2,[1,nframes_track]);
    trx(j).b = reshape(trxb(jj,i,:)/2,[1,nframes_track]);
    trx(j).theta = reshape(trxtheta(jj,i,:),[1,nframes_track]);

    trx(j).firstframe = firstframetrack;
    trx(j).endframe = firstframetrack+nframes_track-1;
    trx(j).nframes = nframes_track;
    trx(j).off = 1-firstframetrack;
    trx(j).roi = i;
    trx(j).roipts = rois{i};
    trx(j).roibb = roibbs(i,:);
    trx(j).moviefile = moviefile;
    trx(j).dt = diff(timestamps(trx(j).firstframe:trx(j).endframe));
    trx(j).timestamp = timestamps(trx(j).firstframe:trx(j).endframe);
  end
  j = sub2ind([2,nrois],1,i);
  trx(j).sex = repmat({'M'},[1,nframes_track]);
  j = sub2ind([2,nrois],2,i);
  trx(j).sex = repmat({'F'},[1,nframes_track]);
end    
outfilename = fullfile(rootdatadir,experiment_name,'trx.mat');
save(outfilename,'trx','timestamps','timestamp_tracked',...
  'moviefile','annfile','bg_nframes',...
  'firstframetrack','lastframetrack',...
  'bgmode','bgthresh','minccarea',...
  'gmmem_params','err_params',...
  'choose_orientations_params','bgmed','rois');

expdir = fullfile(rootdatadir,experiment_name);

%% plot the trajectories
figure(1);
clf;
imagesc(bgmed,[0,255]); axis image;
hold on;
for i = 1:nrois,
  j = i*2-1;
  plot(trx(j).x,trx(j).y,'-','color',[.7,0,0]);
end

title('Trajectories');
savefig(fullfile(expdir,'trajectories.png'),1,'png');

%% plot the heatmap

nbins = 21;

roimuxs = nan(1,nrois);
roimuys = nan(1,nrois);
roiradii = nan(1,nrois);

for roii = 1:nrois,
  tmp = XGRID(roibbs(roii,3):roibbs(roii,4),roibbs(roii,1):roibbs(roii,2));
  tmp(~inrois{roii}) = nan;
  mux = nanmean(tmp(:));
  tmp = YGRID(roibbs(roii,3):roibbs(roii,4),roibbs(roii,1):roibbs(roii,2));
  tmp(~inrois{roii}) = nan;
  muy = nanmean(tmp(:));
  b = bwboundaries(inrois{roii});
  r = sqrt(mean((b{1}(:,2)-mux+roibbs(roii,1)-1).^2+(b{1}(:,1)-muy+roibbs(roii,3)-1).^2));
  roimuxs(roii) = mux;
  roimuys(roii) = muy;
  roiradii(roii) = r;
end

meanradius = mean(roiradii);

binedges = linspace(-meanradius,meanradius,nbins+1);
bincenters = (binedges(1:end-1)+binedges(2:end))/2;

fracim = nan(size(bgmed));

for roii = 1:nrois,
  j = roii*2-1;
  counts = hist3( [trx(j).x(:)-roimuxs(roii),trx(j).y(:)-roimuys(roii)], {bincenters,bincenters} );
  counts = counts';
  frac = counts / sum(counts(:));  
  fracimcurr = max(0,imresize(frac,2*round(meanradius)+[1,1],'nearest'));
  xlims = round(roimuxs(roii))+round(meanradius)*[-1,1];
  ylims = round(roimuys(roii))+round(meanradius)*[-1,1];
  tmp = fracim(ylims(1):ylims(2),xlims(1):xlims(2));
  inroicurr = idxroi(ylims(1):ylims(2),xlims(1):xlims(2))==roii;
  tmp(inroicurr) = fracimcurr(inroicurr);
  fracim(ylims(1):ylims(2),xlims(1):xlims(2)) = tmp;
end

clim = [0,max(fracim(:))];
cminterp = logscale_colormap(jet(1000),[0,clim(2)],(clim(2)-clim(1))/256);
fracim_rgb = colormap_image(fracim,cminterp);
im = repmat(double(bgmed)/255.*(1-double(inroi)/2),[1,1,3])+bsxfun(@times,fracim_rgb,double(inroi)/2);

figure(2);
clf;
him2 = imagesc(0,clim);
hold on;
him = image(im); axis image;
colormap(cminterp);
set(gca,'CLim',clim);
hcb = colorbar;
hfly = nan(1,nrois);
for roii = 1:nrois,
  j = 2*roii;
  hfly(roii) = drawflyo(trx(j),1);
  set(hfly(roii),'Color','w');
end

title('Heatmap of residency time, all frames');
savefig(fullfile(expdir,'heatmap_allframes.png'),2,'png');

%% plot the heatmap, only looking at entrances to new bins

fracim_entrance = nan(size(bgmed));
binsize = mean(diff(binedges));

for roii = 1:nrois,
  j = roii*2-1;
  xidx = floor( (trx(j).x-roimuxs(roii) - binedges(1))/binsize );
  yidx = floor( (trx(j).y-roimuys(roii) - binedges(1))/binsize );
  isentrance = [true,xidx(2:end)~=xidx(1:end-1) | yidx(2:end)~=yidx(1:end-1)];
  
  counts = hist3( [trx(j).x(isentrance)'-roimuxs(roii),trx(j).y(isentrance)'-roimuys(roii)], {bincenters,bincenters} );
  counts = counts';
  frac = counts / sum(counts(:));  
  fracimcurr = max(0,imresize(frac,2*round(meanradius)+[1,1],'nearest'));
  xlims = round(roimuxs(roii))+round(meanradius)*[-1,1];
  ylims = round(roimuys(roii))+round(meanradius)*[-1,1];
  tmp = fracim_entrance(ylims(1):ylims(2),xlims(1):xlims(2));
  inroicurr = idxroi(ylims(1):ylims(2),xlims(1):xlims(2))==roii;
  tmp(inroicurr) = fracimcurr(inroicurr);
  fracim_entrance(ylims(1):ylims(2),xlims(1):xlims(2)) = tmp;
end

clim = [0,max(fracim_entrance(:))];
cminterp = logscale_colormap(jet(1000),clim,(clim(2)-clim(1))/256);
fracim_entrance_rgb = colormap_image(fracim_entrance,cminterp);
im = repmat(double(bgmed)/255.*(1-double(inroi)/2),[1,1,3])+bsxfun(@times,fracim_entrance_rgb,double(inroi)/2);

figure(3);
clf;
him2 = imagesc(0,clim);
hold on;
him = image(im); axis image;
colormap(cminterp);
set(gca,'CLim',clim);
hcb = colorbar;
hfly = nan(1,nrois);
for roii = 1:nrois,
  j = 2*roii;
  hfly(roii) = drawflyo(trx(j),1);
  set(hfly(roii),'Color','w');
end

title('Heatmap of (normalized) number of bin entrances');
savefig(fullfile(expdir,'heatmap_binentrances.png'),3,'png');

%% time series of distance to other fly

roiorder = [1,4,2,5,3,6,7,10,8,11,9,12];
tmp = {'mated','rejected'};
roitype = tmp([1,1,1,2,2,2,2,2,2,1,1,1]);

dcenter = nan(nrois,nframes);
for roii = 1:nrois,
  j1 = 2*roii-1;
  j2 = 2*roii;
  dcenter(roii,:) = sqrt((trx(j1).x-trx(j2).x).^2 + (trx(j1).y-trx(j2).y).^2);  
end
maxdcenter = max(dcenter(:));

figure(4);
clf;
hax = createsubplots(2,ceil(nrois/2),[.05,.01;.1,.05]);
for i = 1:nrois,
  roii = roiorder(i);
  plot(hax(i),dcenter(roii,:),'k.-');
  axis(hax(i),[0,nframes+1,0,maxdcenter*1.05]);
  if mod(i,2) == 0,
    xlabel(hax(i),'Frame');
  else
    set(hax(i),'XTickLabel',{});
  end
  if i > 2,
    set(hax(i),'YTickLabel',{});
  else
    ylabel(hax(i),'Distance apart (centroids)')
  end
  title(hax(i),sprintf('%s (%d)',roitype{roii},roii));
end

savefig(fullfile(expdir,'dcenter_timeseries.png'),4,'png');

%% histogram


nbins = 20;
%binedges = linspace(0,maxdcenter,nbins+1);
%bincenters = (binedges(1:end-1)+binedges(2:end))/2;

[binedges,bincenters] = SelectHistEdges(nbins,[0,maxdcenter],'log');
binsize = diff(binedges);


counts = hist(dcenter',bincenters)';
frac = bsxfun(@rdivide,counts,sum(counts,2));
normfrac = bsxfun(@rdivide,frac,binsize);
maxfrac = max(normfrac(:));

figure(5);
clf;
hax = createsubplots(2,ceil(nrois/2),[.05,.01;.1,.05]);
for i = 1:nrois,
  roii = roiorder(i);
  plot(hax(i),bincenters,normfrac(roii,:),'k.-');
  axis(hax(i),[binedges(1),binedges(end),0,maxfrac*1.05]);
  if mod(i,2) == 0,
    xlabel(hax(i),'Distance apart');
  else
    set(hax(i),'XTickLabel',{});
  end
  if i > 2,
    set(hax(i),'YTickLabel',{});
  else
    ylabel(hax(i),'Fraction of frames / bin size')
  end
  title(hax(i),sprintf('%s (%d)',roitype{roii},roii));
end

savefig(fullfile(expdir,'dcenter_histogram_perroi.png'),5,'png');

%% plot all on one axis

figure(6);
clf;
hold on;
h = nan(1,nrois);
for roii = 1:nrois,
  h(roii) = plot(bincenters,normfrac(roii,:),'.-');
  if strcmpi(roitype{roii},'mated'),
    color = [0,0,0];
  else
    color = [.7,0,0];
  end
  set(h(roii),'Color',color);
end
i = [find(strcmpi(roitype,'mated'),1),find(~strcmpi(roitype,'mated'),1)];
legend(h(i),{'mated','rejected'});
axis([binedges(1),binedges(end),0,maxfrac*1.05]);
xlabel('Distance apart');
ylabel('Fraction of frames / bin size');

meanmated = nanmean(normfrac(strcmpi(roitype,'mated'),:),1);
meanrejected = nanmean(normfrac(~strcmpi(roitype,'mated'),:),1);
stdmated = nanstd(normfrac(strcmpi(roitype,'mated'),:),1,1);
stdrejected = nanstd(normfrac(~strcmpi(roitype,'mated'),:),1,1);
plot(bincenters,meanmated,'-','Color',[0,0,0],'LineWidth',4);
plot(bincenters,meanrejected,'-','Color',[.7,0,0],'LineWidth',4);

savefig(fullfile(expdir,'dcenter_histogram.png'),6,'png');
