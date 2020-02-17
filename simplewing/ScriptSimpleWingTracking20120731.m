%% set up paths

realrootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
rootdatadir = '/groups/branson/bransonlab/projects/olympiad/WingTrackedData';
addpath ../misc;
addpath ../filehandling;
addpath /groups/branson/home/bransonk/olympiad_bransonlab/FlyBowlAnalysis;

%% parameters

experiment_names = {
  'GMR_53E04_AE_01_TrpA_Rig1Plate15BowlA_20120201T113036'
  'GMR_53E04_AE_01_TrpA_Rig1Plate15BowlB_20120201T113031'
  'GMR_53E04_AE_01_TrpA_Rig1Plate15BowlC_20120201T113025'
  'GMR_53E04_AE_01_TrpA_Rig1Plate15BowlD_20120201T113021'
  'GMR_10H09_AE_01_TrpA_Rig1Plate10BowlA_20110728T112504'
  'GMR_10H09_AE_01_TrpA_Rig1Plate10BowlB_20110728T112510'
  'GMR_10H09_AE_01_TrpA_Rig1Plate10BowlC_20110728T112327'
  'GMR_10H09_AE_01_TrpA_Rig1Plate10BowlD_20110728T112331'
  'GMR_46E03_AE_01_TrpA_Rig1Plate10BowlA_20110624T144217'
  'GMR_46E03_AE_01_TrpA_Rig1Plate10BowlB_20110624T144212'
  'GMR_46E03_AE_01_TrpA_Rig1Plate10BowlC_20110624T144055'
  'GMR_46E03_AE_01_TrpA_Rig1Plate10BowlD_20110624T144052'
  };

%experiment_name = 'pBDPGAL4U_TrpA_Rig2Plate17BowlA_20110929T143440';
moviefilestr = 'movie.ufmf';
annfilestr = 'movie.ufmf.ann';
trxfilestr = 'registered_trx.mat';
outtrxfilestr = 'wingtracking_results.mat';
DEBUG = false;

mindwing_high = 50;
mindwing_low = 30;
mindbody = 100;
radius_dilate_body = 1;
radius_open_wing = 1;
min_wingcc_area = 3;
max_wingcc_dist = .25;
max_wingpx_angle = 3*pi/4;
max_wing_angle_otherside = pi/18;
min_nonzero_wing_angle = pi/18;
wing_min_peak_dist_bins = 3;
wing_min_peak_threshold_frac = 0;
nbins_dthetawing = 50;
edges_dthetawing = linspace(-max_wingpx_angle,max_wingpx_angle,nbins_dthetawing+1);
centers_dthetawing = (edges_dthetawing(1:end-1)+edges_dthetawing(2:end))/2;
binsize_dthetawing = mean(diff(centers_dthetawing));
wing_peak_min_frac = 1/nbins_dthetawing*2;
wing_frac_filter = [1,2,1];
wing_frac_filter = wing_frac_filter / sum(wing_frac_filter);
wing_fit_method = 'peaks';
wing_min_prior = .25;
min_single_wing_area = 10;
firstframe = 1;

se_dilate_body = strel('disk',radius_dilate_body);
se_open_wing = strel('disk',radius_open_wing);

% for sub-bin accuracy in fitting the wings
wing_radius_quadfit_bins = 1;
subbin_x = (-wing_radius_quadfit_bins:wing_radius_quadfit_bins)';    

units = struct(...
  'nwingsdetected',parseunits('unit'),...
  'wing_areal',parseunits('px^2'),...
  'wing_arear',parseunits('px^2'),...
  'wing_trough_angle',parseunits('rad'));

%% copy over data

for expi = 1:numel(experiment_names),
  expdir = fullfile(realrootdatadir,experiment_names{expi});
  SymbolicCopyExperimentDirectory(expdir,rootdatadir);  
end

%% loop over experiments

for expi = 1:numel(experiment_names),
  
experiment_name = experiment_names{expi};

%% read stuff
expdir = fullfile(rootdatadir,experiment_name);

moviefile = fullfile(expdir,moviefilestr);
annfile = fullfile(expdir,annfilestr);
trxfile = fullfile(expdir,trxfilestr);
[readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefile);
bgmodel = read_ann(annfile,'background_center');
bgmodel = reshape(bgmodel,[headerinfo.nc,headerinfo.nr])';
[trx,~,~,timestamps] = load_tracks(trxfile);

nr = headerinfo.nr;
nc = headerinfo.nc;
npx = headerinfo.nr*headerinfo.nc;

[XGRID,YGRID] = meshgrid(1:nc,1:nr);

clear perframedata;

%% start tracking

if DEBUG,
  colors = hsv(numel(trx));
  colors = colors(randperm(numel(trx)),:);
  hims = nan(1,4);
  if DEBUG > 1,
    hfig = 1;
    figure(hfig);
    clf;
    hax = createsubplots(1,3,.01);
  end
  drawnow;
end

%for t = t:max([trx.endframe]),
for t = min([trx.firstframe]):max([trx.endframe]),
%for t = round(linspace(min([trx.firstframe]),max([trx.endframe]),20)),

if mod(t,30) == 0,
  fprintf('Frame %d\n',t);
  drawnow;
end
  
im = double(readframe(t));
dbkgd = bgmodel - im;

%% morphology to get foreground, body wing pixels

isbody_thresh = dbkgd >= mindbody;
isbody = imdilate(isbody_thresh,se_dilate_body);

% isfore_low = dbkgd >= mindwing_low;
% isfore_high = dbkgd >= mindwing_high;

iswing_high = ~isbody & dbkgd >= mindwing_high;
iswing_low = ~isbody & dbkgd >= mindwing_low;
iswing = imreconstruct(iswing_high,iswing_low,4);
iswing = imopen(iswing,se_open_wing);
iswing = imclose(iswing,se_open_wing);

isfore_thresh = imclose(isbody_thresh | iswing,se_dilate_body);

%isfore_thresh = dbkgd > mindwing_high;
idxfore_thresh = find(isfore_thresh);
npxfore_thresh = nnz(isfore_thresh);

fore2dbkgd = dbkgd(isfore_thresh);
fore2body = fore2dbkgd >= mindbody;

body2fore = find(fore2body);
npxbody_thresh = nnz(isbody_thresh);


%% plot

if DEBUG > 1,
  axcurr = 1;
  if isnan(hims(axcurr)),
    cla(hax(axcurr));
  end
  imtmp = repmat(im(:),[1,3]);
  imtmp(isbody,1) = min(imtmp(isbody,1)+100,255);
  imtmp(iswing,2) = min(imtmp(iswing,2)+100,255);
  imtmp = uint8(reshape(imtmp,[nr,nc,3]));
  if isnan(hims(1)),
    hims(axcurr) = image(imtmp,'Parent',hax(axcurr)); axis(hax(axcurr),'image','off');
  else
    set(hims(axcurr),'CData',imtmp);
  end    
  drawnow;
  title(hax(axcurr),sprintf('Segmentation of frame %d into bg, body and wing',t));
end

%% assign targets to connected components
[L,ncc] = bwlabel(isfore_thresh);
fore2cc = L(isfore_thresh);
fly2cc = nan(1,numel(trx));
pcurr = [];
for fly = 1:numel(trx),
  
  if t < trx(fly).firstframe || t > trx(fly).endframe,
    continue;
  end
  i = trx(fly).off+t; 
  x = round(trx(fly).x(i));
  y = round(trx(fly).y(i));
  if y >= 1 && y <= nr && x >= 1 && x <= nc && ...
      isfore_thresh(y,x),
    fly2cc(fly) = L(y,x);
  else
    if isempty(pcurr),
      pcurr = [XGRID(isfore_thresh),YGRID(isfore_thresh)];
    end
    mu = [trx(fly).x(i),trx(fly).y(i)];
    S = axes2cov(trx(fly).a(i)*2,trx(fly).b(i)*2,trx(fly).theta(i));
    diffs = bsxfun(@minus,pcurr,mu);
    c = chol(S);
    temp = diffs/c;
    dcurr = sum(temp.^2, 2);
    [~,j] = min(dcurr);
    fly2cc(fly) = fore2cc(j);
  end
  
end

%% assign unassigned connected components to flies

fly2nccs = ones(1,numel(trx));
unassignedcc = setdiff(1:ncc,fly2cc);
if ~isempty(unassignedcc),
  mus = nan(numel(trx),2);
  for fly = 1:numel(trx),    
    i = trx(fly).off+t;
    if t < trx(fly).firstframe || t > trx(fly).endframe,
      continue;
    end
    mus(fly,:) = [trx(fly).x(i),trx(fly).y(i)];
  end
  for cci = unassignedcc,
    % remove small connected components
    if nnz(L==cci) < min_wingcc_area,
      continue;
    end
    x = XGRID(L==cci);
    y = YGRID(L==cci);
    D = dist2([x(:),y(:)],mus);
    d = min(D,[],1);
    [mind,fly] = min(d);
    i = trx(fly).off+t;
    mind = mind/(trx(fly).a(i)*4);
    if mind > max_wingcc_dist,
      continue;
    end
    L(L==cci) = fly2cc(fly);
    fore2cc(fore2cc==cci) = fly2cc(fly);
    fly2nccs(fly) = fly2nccs(fly) + 1;
  end
end

%% plot

if DEBUG > 1,

axcurr = 2;
if isnan(hims(axcurr)),
  cla(hax(axcurr));
end

imtmp = repmat(im(:),[1,3]);
for fly = 1:numel(trx),
  idx1 = L == fly2cc(fly);
  imtmp(idx1,:) = min(bsxfun(@plus,imtmp(idx1,:)*3,255*colors(fly,:))/4,255);
end
imtmp = uint8(reshape(imtmp,[nr,nc,3]));

if isnan(hims(axcurr)),
  hims(axcurr) = image(imtmp,'Parent',hax(axcurr)); axis(hax(axcurr),'image','off');
else
  set(hims(axcurr),'CData',imtmp);
  delete(htext{axcurr}(ishandle(htext{axcurr})));
end
htext{axcurr} = [];

for cci = unique(fly2cc),
  flies = find(fly2cc == cci);
  s = sprintf('%d,',flies);
  s = s(1:end-1);
  x = mean(XGRID(L==cci));
  y = mean(YGRID(L==cci));
  htext{axcurr} = [htext{axcurr},text(x,y,s,'HorizontalAlignment','center','VerticalAlignment','middle','Clipping','on','Parent',hax(axcurr))];
end

title(hax(axcurr),sprintf('Connected components for frame %d',t));
drawnow;
end

%% loop over connected components, assign pixels to flies

unique_cc = unique(fly2cc(~isnan(fly2cc)));
% fore2fly corresponds to isfore_thresh
fore2fly = nan(npxfore_thresh,1);

xgrid_isfore = XGRID(isfore_thresh);
ygrid_isfore = YGRID(isfore_thresh);

% for watershed
hy = fspecial('sobel');
hx = hy';
se_boundary = strel('disk',1);

for cci = unique_cc,
  flies = find(fly2cc == cci);
  fore2curr = fore2cc == cci;
  % all pixels belong to one fly
  if numel(flies) == 1,
    fore2fly(fore2curr) = flies;
    continue;    
  end
  
  % assign body pixels based on mahdist
  
  % indices of foreground pixels in the current cc
  curr2fore = find(fore2curr);
  ncurr = numel(curr2fore);
  
  % ncurr x 1, whether current cc pixels are body or not
  curr2body = fore2body(fore2curr);
  
  % nbodycurr x 1, indices of current cc pixels in body
  body2curr = find(curr2body);
  nbodycurr = numel(body2curr);
  
  % compute distance from each current body pixel to fly
  mind = inf(nbodycurr,1);
  body2id = nan(nbodycurr,1);
  pbody = [xgrid_isfore(curr2fore(body2curr)),ygrid_isfore(curr2fore(body2curr))];
  for fly = flies,
    i = trx(fly).off+t;
    mu = [trx(fly).x(i),trx(fly).y(i)];
    S = axes2cov(trx(fly).a(i)*2,trx(fly).b(i)*2,trx(fly).theta(i));
    diffs = bsxfun(@minus,pbody,mu);
    c = chol(S);
    temp = diffs/c;
    dcurr = sum(temp.^2, 2);

    idx = dcurr < mind;
    body2id(idx) = fly;
    mind(idx) = dcurr(idx);
  end
  
  % compute distance from each non-body pixel to each body pixel
  pcurr = [xgrid_isfore(fore2curr),ygrid_isfore(fore2curr)];
  xlims = [floor(min(pcurr(:,1)))-1,ceil(max(pcurr(:,1)))+1];
  ylims = [floor(min(pcurr(:,2)))-1,ceil(max(pcurr(:,2)))+1];
  xlims = max(1,min(nc,xlims));
  ylims = max(1,min(nr,ylims));
  dx = diff(xlims)+1;
  dy = diff(ylims)+1;
  
  % pixels we can go through
  bw = L(ylims(1):ylims(2),xlims(1):xlims(2)) == cci;
  npxcurr = nnz(bw);
  mind = inf(size(bw));
  flycurr = nan(size(bw));
  for fly = flies,
    % pixels we are finding the distance to
    c = pbody(body2id==fly,1)-xlims(1)+1;
    r = pbody(body2id==fly,2)-ylims(1)+1;
    
    % find distance within foreground pixels
    dcurr = bwdistgeodesic(bw,c,r,'quasi-euclidean');
    
    % assign if smaller
    tmpidx = dcurr < mind;
    mind(tmpidx) = dcurr(tmpidx);
    flycurr(tmpidx) = fly;
  end
  
  % superpixel segmentation
  imbb = im(ylims(1):ylims(2),xlims(1):xlims(2));
  imbb(~isfore_thresh(ylims(1):ylims(2),xlims(1):xlims(2))) = ...
    mean(imbb(~isfore_thresh(ylims(1):ylims(2),xlims(1):xlims(2))));
  Iy = imfilter(imbb, hy, 'replicate');
  Ix = imfilter(imbb, hx, 'replicate');
  gradmag = sqrt(Ix.^2 + Iy.^2);
  isboundary = ~bw & imdilate(bw,se_boundary);
  gradmag(isboundary) = inf;
  
  Lseg = watershed(gradmag);
  
  % remove bg segments
  nsegcurr = max(Lseg(:));
  for segi = 1:nsegcurr,
    segcurr = Lseg==segi;
    tmp = bw(segcurr);
    tmp(isnan(tmp)) = [];
    if isempty(tmp),
      continue;
    end
    if nnz(tmp) < numel(tmp)/2,
      Lseg(segcurr) = nan;
    end
  end
  
  imbb(~bw|isnan(Lseg)) = nan;

  % assign watershed pixels based on intensity difference
  [rw,cw] = find(Lseg==0);
  iw = sub2ind(size(Lseg),rw,cw);
  nwatershed = numel(rw);
  mind = inf(nwatershed,1);
  Lw = nan(nwatershed,1);
  for dc = -1:1,
    for dr = -1:1,
      if dc == 0 && dr == 0,
        continue;
      end
      % adjacent pixel
      r1 = rw + dr;
      c1 = cw + dc;
      % make sure it is in bounds
      idxw = find(r1 >= 1 & r1 <= dy & c1 >= 1 & c1 <= dx);
      iw1 = sub2ind(size(Lseg),r1(idxw),c1(idxw));
      % also make sure it does not have Lseg == 0, and is foreground
      idxw(Lseg(iw1)==0 | ~bw(iw1)) = [];
      iw1 = sub2ind(size(Lseg),r1(idxw),c1(idxw));
      % compute intensity difference
      dcurr = nan(nwatershed,1);
      lcurr = nan(nwatershed,1);
      dcurr(idxw) = abs(imbb(iw(idxw))-imbb(iw1));
      lcurr(idxw) = Lseg(iw1);
      idxw1 = dcurr < mind;
      if any(isnan(lcurr(idxw1))),
        keyboard;
      end
      %fprintf('dr = %d, dc = %d, assigning for: %s\n',dr,dc,mat2str(find(idxw1)));
      mind(idxw1) = dcurr(idxw1);
      Lw(idxw1) = lcurr(idxw1);
    end
  end
  
  Lseg(iw) = Lw;
  Lseg(~bw) = nan;
  
  % for each superpixel, take the majority vote  
  nsegcurr = max(Lseg(:));
  flycurr1 = nan(size(flycurr));
  for segi = 1:nsegcurr,
    segcurr = Lseg==segi;
    tmp = flycurr(segcurr);
    tmp(isnan(tmp)) = [];
    if isempty(tmp),
      continue;
    end
    flyseg = mode(tmp);
    flycurr1(segcurr) = flyseg;
  end
  
  flycurr1(~bw|isnan(flycurr1)) = 0;
  
  % store
  fore2fly(fore2curr) = flycurr1(bw);
  
end

%% plot segmentation of one connected component into flies
% 
% tmpplot = repmat(imbb(:),[1,3]);
% colors = lines(numel(flies));
% colors = colors(randperm(size(colors,1)),:);
% for i = 1:numel(flies),
%   idx1 = flycurr1==flies(i);
%   tmpplot(idx1,:) = min(bsxfun(@plus,tmpplot(idx1,:),255*colors(i,:))/2,255);
% end
% 
% image(uint8(reshape(tmpplot,[size(imbb),3]))); axis(hax(axcurr),'image','off')

%% plot segmentation of one connected component into superpixels
% 
% tmpplot = repmat(imbb(:),[1,3]);
% nseg = double(max(Lseg(:)))+1;
% colors = hsv(nseg);
% colors = colors(randperm(size(colors,1)),:);
% for i = 1:nseg+1,
%   idx1 = Lseg==i-1;
%   if ~any(idx1(:)),
%     continue;
%   end
%   tmpplot(idx1,:) = min(bsxfun(@plus,tmpplot(idx1,:)*3,255*colors(i,:))/4,255);
% end
% 
% clf;
% image(uint8(reshape(tmpplot,[size(imbb),3]))); axis(hax(axcurr),'image','off')
% hold on;
% for i = 1:nseg+1,
%   [rtmp,ctmp] = find(Lseg==i-1);
%   if ~isempty(rtmp),
%     text(mean(ctmp),mean(rtmp),num2str(i-1),'horizontalalignment','center');
%   end
% end

%% plot

if DEBUG > 1,
  
axcurr = 3;
if isnan(hims(axcurr)),
  cla(hax(axcurr));
end

imtmp = repmat(im(:),[1,3]);
for fly = 1:numel(trx),
  idx1 = idxfore_thresh(fore2fly==fly);
  imtmp(idx1,:) = min(bsxfun(@plus,imtmp(idx1,:)*3,255*colors(fly,:))/4,255);
end

imtmp = uint8(reshape(imtmp,[nr,nc,3]));
if isnan(hims(axcurr)),
  hims(axcurr) = image(imtmp,'Parent',hax(axcurr)); axis(hax(axcurr),'image','off');
else
  set(hims(axcurr),'CData',imtmp);
  delete(htext{axcurr}(ishandle(htext{axcurr})));
end
htext{axcurr} = [];

for fly = 1:numel(trx),
  s = sprintf('%d',fly);
  idx1 = idxfore_thresh(fore2fly==fly);
  x = mean(XGRID(idx1));
  y = mean(YGRID(idx1));
  htext{axcurr} = [htext{axcurr},text(x,y,s,'HorizontalAlignment','center','VerticalAlignment','middle','Clipping','on','Parent',hax(axcurr))];
end

title(hax(axcurr),sprintf('Segmentation for frame %d',t));
drawnow;

end

%% find pixels that belong to each fly's wings
% then fit wings

if ~exist('perframedata','var'),
  perframedata = struct;
end
if ~isfield(perframedata,'nwingsdetected'),
  perframedata.nwingsdetected = cell(1,numel(trx));
end
if ~isfield(perframedata,'wing_areal'),
  perframedata.wing_areal = cell(1,numel(trx));
end
if ~isfield(perframedata,'wing_arear'),
  perframedata.wing_arear = cell(1,numel(trx));
end
if ~isfield(perframedata,'wing_trough_angle'),
  perframedata.wing_trough_angle = cell(1,numel(trx));
end
if ~isfield(trx,'wing_anglel'),
  for fly = 1:numel(trx),
    trx(fly).wing_anglel = [];
  end
end
if ~isfield(trx,'wing_angler'),
  for fly = 1:numel(trx),
    trx(fly).wing_angler = [];
  end
end

fore2wing = iswing(isfore_thresh);
fore2flywing = zeros(size(fore2fly));
for fly = 1:numel(trx),

  if t < trx(fly).firstframe || t > trx(fly).endframe,
    continue;
  end
  i = trx(fly).off+t; 
  x = trx(fly).x(i);
  y = trx(fly).y(i);
  theta = trx(fly).theta(i);
  idxcurr = fore2wing&fore2fly==fly;
  xwing = xgrid_isfore(idxcurr);
  ywing = ygrid_isfore(idxcurr);
  dthetawing = modrange(atan2(ywing-y,xwing-x)-(theta+pi),-pi,pi);
  isallowed = abs(dthetawing) <= max_wingpx_angle;
  idxcurr = find(idxcurr);
  idxcurr(~isallowed) = [];
  fore2flywing(idxcurr) = fly;
  dthetawing(~isallowed) = [];
  nwingpx = numel(dthetawing);
  
  % fit wings
  
  switch wing_fit_method,
    
    case 'peaks',

      frac = histc(dthetawing,edges_dthetawing)/nwingpx;
      frac = frac(:)';
      frac = [frac(1:end-2),frac(end-1)+frac(end)];
      smoothfrac = conv(frac,wing_frac_filter,'same');
      %smoothfrac = imfilter(frac,wing_frac_filter);      
      % wings are at local maxima
      
      % this is much faster than calling the findpeaks function
      locs = []; pks = [];
      if nwingpx > min_single_wing_area,
        [pk,loc] = max(smoothfrac);
        if pk >= wing_min_peak_threshold_frac,
          locs(1) = loc; pks(1) = pk;
          isallowed = [true,smoothfrac(2:end)>smoothfrac(1:end-1)] & ...
            [smoothfrac(1:end-1)>=smoothfrac(2:end),true];
          isallowed(max(1,loc-wing_min_peak_dist_bins):min(nbins_dthetawing,loc+wing_min_peak_dist_bins)) = false;
          [pk,loc] = max(smoothfrac(isallowed));
          if pk >= wing_peak_min_frac,
            idxallowed = find(isallowed);
            locs(2) = idxallowed(loc);
            pks(2) = pk;
          end
        end
      end
            
%       [pks,locs] = findpeaks(smoothfrac,'MinPeakHeight',wing_peak_min_frac,...
%         'MinPeakDistance',wing_min_peak_dist_bins,...
%         'Threshold',wing_min_peak_threshold_frac,...
%         'NPeaks',2,...
%         'SortStr','descend');
      
      % sub-bin precision
      if isempty(locs),
        wing_angles_curr = 0;
        npeaks = 0;
        area = 0;
        wing_trough_angle = 0;
      elseif numel(locs) == 1,
         wing_angles_curr = median(dthetawing);
         npeaks = 1;
         area = nwingpx;
         wing_trough_angle = wing_angles_curr;
      else
        wing_angles_curr = centers_dthetawing(locs);
        npeaks = 2;
        
        for j = 1:numel(locs),
          subbin_x_curr = subbin_x+locs(j);
          subbin_x_curr(subbin_x_curr < 1 | subbin_x_curr > nbins_dthetawing) = [];
          ncurr = numel(subbin_x_curr);
          Xcurr = [ones(ncurr,1),  subbin_x_curr,  subbin_x_curr.^2];
          
          % frac = X*coeffs
          coeffs = Xcurr\smoothfrac(subbin_x_curr)';
          
          % max(c1 + c2*x + c3*x^2)
          % c2 + 2*c3*x = 0
          % x = -c2 / (2*c3)
          
          % make sure it is concave
          if coeffs(3) >= 0,
            continue;
          end
          
          maxx = -coeffs(2)/(2*coeffs(3));
          
          % make sure it is within range
          if abs(maxx-locs(j)) > 1 || maxx > nbins_dthetawing || maxx < 1,
            continue;
          end
          
          % interpolate
          wceil = maxx-floor(maxx);
          wing_angles_curr(j) = centers_dthetawing(floor(maxx))*(1-wceil) + wceil*centers_dthetawing(ceil(maxx));
          
        end
      
        % can't have two wings on the same side of the body
        if sign(wing_angles_curr(1)) == sign(wing_angles_curr(2)) && ...
            min(abs(wing_angles_curr)) >= min_nonzero_wing_angle,
          wing_angles_curr(end) = [];
          npeaks = 1;
          area = nwingpx;
        else
          
          % find the trough between the two peaks
          minloc = min(locs);
          maxloc = max(locs);
          [~,troughloc] = min(smoothfrac(minloc:maxloc));
          troughloc = troughloc + minloc - 1;
          % find other locations with the same value
          loc1 = find(smoothfrac(minloc:troughloc)>smoothfrac(troughloc),1,'last');
          if isempty(loc1),
            loc1 = troughloc;
          else
            loc1 = loc1+minloc;
          end
          loc2 = find(smoothfrac(troughloc:maxloc)>smoothfrac(troughloc),1,'first');
          if isempty(loc2),
            loc2 = troughloc;
          else
            loc2 = loc2+troughloc-2;
          end
          troughloc = (loc1+loc2)/2;
          if mod(troughloc,1) > 0,
            area = sum(smoothfrac(1:troughloc-.5));
            wing_trough_angle = (centers_dthetawing(troughloc-.5)+centers_dthetawing(troughloc+.5))/2;
          else
            area = sum(smoothfrac(1:troughloc-1))+smoothfrac(troughloc)/2;
            wing_trough_angle = centers_dthetawing(troughloc);
          end
          area = [area,1-area]*nwingpx;
          
          if all(area < min_single_wing_area),
            wing_angles_curr = 0;
            npeaks = 0;
            area = 0;
            wing_trough_angle = 0;
          elseif area(1) < min_single_wing_area,
            [~,removei] = min(wing_angles_curr);
            area(1) = [];
            wing_angles_curr(removei) = [];
            npeaks = 1;
          elseif area(2) < min_single_wing_area,
            [~,removei] = max(wing_angles_curr);
            area(2) = [];
            wing_angles_curr(removei) = [];
            npeaks = 1;
          end
        end
        
      end
      wing_angles_curr = sort(wing_angles_curr);
      if numel(wing_angles_curr) == 1,
        if abs(wing_angles_curr) > min_nonzero_wing_angle,
          wing_trough_angle = 0;
          [wing_angles_curr,order] = sort([0,wing_angles_curr]);
          area = [0,area]; %#ok<AGROW>
          area = area(order);
        else
          wing_angles_curr = repmat(wing_angles_curr,[1,2]);
          area = repmat(area,[1,2]);
        end
      end
      perframedata.nwingsdetected{fly}(i) = npeaks;
      perframedata.wing_areal{fly}(i) = area(1);
      perframedata.wing_arear{fly}(i) = area(2);
      perframedata.wing_trough_angle{fly}(i) = wing_trough_angle;
      
    case 'gmm',
  
      gmmweights = max(0,fore2dbkgd(idxcurr)-mindwing_high);
      
      % or try gmm-based approach
      aic = nan(1,2);
      bic = nan(1,2);
      mu = cell(1,2);
      S = cell(1,2);
      priors = cell(1,2);
      
      k = 1;
      [mu{k},S{k}] = weighted_mean_cov(dthetawing,gmmweights);
      priors{k} = 1;
      nllcurr = -sum(log(normpdf(dthetawing,mu{k},sqrt(S{k}))).*gmmweights);
      aic(k) = 2*nllcurr+2*3*k;
      bic(k) = 2*nllcurr+3*k*log(nwingpx);
      
      k = 2;
      if t == firstframe || i == 1,
        [mucurr,Scurr,priorscurr,postcurr,nllcurr] = mygmm(dthetawing,k,'Replicates',10,'weights',gmmweights);
      else
        startmu = [trx(fly).wing_anglel(i-1);trx(fly).wing_angler(i-1)];
        dmu = abs(startmu(1)-startmu(2));
        if dmu < .02;
          startmu(1) = startmu(1) - (.01-dmu);
          startmu(2) = startmu(2) + (.01-dmu);
        end
        [mucurr,Scurr,priorscurr,postcurr,nllcurr] = mygmm(dthetawing,k,'Start',startmu,'weights',gmmweights);
      end
      mu{k} = mucurr; S{k} = Scurr; priors{k} = priorscurr;
      aic(k) = 2*nllcurr+2*3*k;
      bic(k) = 2*nllcurr+3*k*log(sum(gmmweights));
%       if min(priors{2}) < wing_min_prior,
%         k = 1;
%       else
%         [~,k] = min(bic);
%       end
      [~,k] = min(aic);
      wing_angles_curr = mu{k};
      if numel(wing_angles_curr) == 1,
        wing_angles_curr = repmat(wing_angles_curr,[1,2]);
      end
  end
  
  trx(fly).wing_anglel(i) = wing_angles_curr(1);
  trx(fly).wing_angler(i) = wing_angles_curr(2);
  
end

%% plot

if DEBUG,
  
  hfig = 2;
  axcurr = 4;
  figure(hfig);
  clf;
  hax1 = gca;
  imtmp = repmat(im(:),[1,3]);
  for fly = 1:numel(trx),
    idx1 = idxfore_thresh(fore2flywing==fly);
    imtmp(idx1,:) = min(bsxfun(@plus,imtmp(idx1,:)*3,255*colors(fly,:))/4,255);
  end
  
  if isnan(hims(axcurr)) || ~ishandle(hims(axcurr)),
    cla(hax1);
    hims(axcurr) = image(uint8(reshape(imtmp,[nr,nc,3])),'Parent',hax1);
    axis(hax1,'image','off');
    hold(hax1,'on');
    if DEBUG > 1,
      linkaxes([hax,hax1]);
    end
  else
    set(hims(axcurr),'CData',uint8(reshape(imtmp,[nr,nc,3])));
  end
  
  if exist('hwing','var'),
    delete(hwing(ishandle(hwing)));
  end
  if exist('htext2','var'),
    delete(htext2(ishandle(htext2)));
  end
  if exist('htrough','var'),
    delete(htrough(ishandle(htrough)));
  end
  hwing = [];
  htrough = [];
  htext2 = [];
  for fly = 1:numel(trx),
    if t < trx(fly).firstframe || t > trx(fly).endframe,
      continue;
    end
    i = trx(fly).off+t;
    xwing = [nan,trx(fly).x(i),nan];
    ywing = [nan,trx(fly).y(i),nan];
    wing_angles = [trx(fly).wing_anglel(i),trx(fly).wing_angler(i)];
    xwing([1,3]) = trx(fly).x(i) + 4*trx(fly).a(i)*cos(trx(fly).theta(i)+pi+wing_angles);
    ywing([1,3]) = trx(fly).y(i) + 4*trx(fly).a(i)*sin(trx(fly).theta(i)+pi+wing_angles);
    xtrough = trx(fly).x(i)+2*trx(fly).a(i)*cos(trx(fly).theta(i)+pi+perframedata.wing_trough_angle{fly}(i));
    ytrough = trx(fly).y(i)+2*trx(fly).a(i)*sin(trx(fly).theta(i)+pi+perframedata.wing_trough_angle{fly}(i));
    hwing(end+1) = plot(hax1,xwing,ywing,'.-','color',colors(fly,:));
    htrough(end+1) = plot(hax1,xtrough,ytrough,'x','color',colors(fly,:));
    htext2(end+1) = text(xwing(1),ywing(1),sprintf('%.1f',perframedata.wing_areal{fly}(i)));
    htext2(end+1) = text(xwing(3),ywing(3),sprintf('%.1f',perframedata.wing_arear{fly}(i)));
  end
  
  if exist('htext','var') && numel(htext) >= axcurr,
    delete(htext{axcurr}(ishandle(htext{axcurr})));
  end
  
  htext{axcurr} = [];
  
  for fly = 1:numel(trx),
    s = sprintf('%d',fly);
    idx1 = idxfore_thresh(fore2fly==fly);
    x = mean(XGRID(idx1));
    y = mean(YGRID(idx1));
    htext{axcurr} = [htext{axcurr},text(x,y,s,'HorizontalAlignment','center','VerticalAlignment','middle','Clipping','on','Parent',hax1,'Color','w')];
  end
  
  if DEBUG > 1,    
    input('');
%     tmpoutfilename = sprintf('WingTrackingExamples_Segmentation_%05d.pdf',t);
%     savefig(tmpoutfilename,1,'pdf');
%     tmpoutfilename = sprintf('WingTrackingExamples_Results_%05d.pdf',t);
%     savefig(tmpoutfilename,2,'pdf');
  else
    drawnow;
  end
end

if mod(t,5000) == 0,
  save tmp.mat trx perframedata t
end


end

%%

outtrxfile = fullfile(expdir,outtrxfilestr);
save_tracks(trx,outtrxfile,'timestamps',timestamps);

%%
% trxfile = fullfile(expdir,outtrxfilestr);
% load(trxfile);

for fly = 1:numel(trx),    
  trx(fly).xwingl = trx(fly).x + 4*trx(fly).a.*cos(trx(fly).theta+ pi+trx(fly).wing_anglel);
  trx(fly).ywingl = trx(fly).y + 4*trx(fly).a.*sin(trx(fly).theta+ pi+trx(fly).wing_anglel);
  trx(fly).xwingr = trx(fly).x + 4*trx(fly).a.*cos(trx(fly).theta+ pi+trx(fly).wing_angler);
  trx(fly).ywingr = trx(fly).y + 4*trx(fly).a.*sin(trx(fly).theta+ pi+trx(fly).wing_angler);
end

save(outtrxfile,'trx','timestamps');

fns = fieldnames(perframedata);
for i = 1:numel(fns),
  fn = fns{i};
  s = struct('data',{perframedata.(fn)},'units',units.(fn));
  filename = fullfile(expdir,'perframe',[fn,'.mat']);
  save(filename,'-struct','s');
end


%%

try
  fclose(fid);
end

end