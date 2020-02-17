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
rootdatadir = '/groups/branson/bransonlab/projects/JAABA/data/eric_courtshipbowls';
experiment_name = 'CantonS_TrpA_Rig2Plate17BowlC_20121207T142430';
%rootdatadir = 'C:\Data\CourtshipBowl';
%experiment_name = '20120719_shelby';

addpath /groups/branson/bransonlab/projects/CourtshipBowls/CourtshipBowlAnalysis
settingsdir = '/groups/branson/bransonlab/projects/CourtshipBowls/CourtshipBowlAnalysis/settings';
analysis_protocol = 'current';


%% parameters

moviefile = fullfile(rootdatadir,experiment_name,'movie.ufmf');
annfile = fullfile(rootdatadir,experiment_name,'movie.ufmf.ann');

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
bgthresh = 120;

% minimum area of a connected component in pixels
minccarea = 20;

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
  buffer = repmat(buffer,[1,1,bg_nframes]);
  frames = round(linspace(1,nframes,bg_nframes));
  for i = 1:bg_nframes,
    t = frames(i);
    buffer(:,:,i) = readframe(t);
  end
  
  bgmed = uint8(median(single(buffer),3));
end

%% get regions of interest

roiparams = ReadParams(fullfile(settingsdir,analysis_protocol,'roidetect_params.txt'));
nrois = floor(numel(roiparams.roimus)/2);
roiparams.roimus = reshape(roiparams.roimus(1:2*nrois),[2,nrois])';
roidata = struct;
[roidata.centerx,roidata.centery,roidata.radii,roidata.detectscores] = ...
  detectCourtshipBowlROIs(bgmed,roiparams);

% create masks
[XGRID,YGRID] = meshgrid(1:imwidth,1:imheight);
roibbs = [max(1,floor(roidata.centerx(:)-roidata.radii(:))),...
  min(imwidth,ceil(roidata.centerx(:)+roidata.radii(:))),...
  max(1,floor(roidata.centery(:)-roidata.radii(:))),...
  min(imheight,ceil(roidata.centery(:)+roidata.radii(:)))];

idxroi = zeros(imheight,imwidth);
inrois = cell(1,nrois);
for i = 1:nrois,
  bb = roibbs(i,:);
  inrois{i} = (XGRID(bb(3):bb(4),bb(1):bb(2)) - roidata.centerx(i)).^2 + ...
    (YGRID(bb(3):bb(4),bb(1):bb(2)) - roidata.centery(i)).^2 ...
    <= roidata.radii(i)^2;
  tmp = idxroi(bb(3):bb(4),bb(1):bb(2));
  tmp(inrois{i}) = i;
  idxroi(bb(3):bb(4),bb(1):bb(2)) = tmp;
end
roidata.roibbs = roibbs;
roidata.idxroi = idxroi;
roidata.inrois = inrois;

% do background subtraction to count flies in each roi
framessample = round(linspace(1,nframes,roiparams.nframessample));
areassample = cell(nrois,roiparams.nframessample);

for i = 1:roiparams.nframessample,
  im = readframe(framessample(i));
  switch bgmode,
    case DARKBKGD,
      dbkgd = imsubtract(im,bgmed);
    case LIGHTBKGD,
      dbkgd = imsubtract(bgmed,im);
    case OTHERBKGD,
      dbkgd = imabsdiff(im,bgmed);
  end
  
  % threshold
  isfore = dbkgd >= bgthresh;
  
  for j = 1:nrois,
    roibb = roibbs(j,:);
    isforebb = isfore(roibb(3):roibb(4),roibb(1):roibb(2));
    isforebb(~inrois{j}) = false;
    cc = bwconncomp(isforebb);
    areassample{j,i} = cellfun(@numel,cc.PixelIdxList);    
  end
  
end

nflies_per_roi = nan(1,nrois);
for i = 1:nrois,
  nccs = cellfun(@(x) nnz(x >= minccarea),areassample(i,:));
  max_nccs = max(nccs);
  mode_nccs = mode(nccs);
  if mode_nccs ~= 1,
    nflies_per_roi(i) = mode_nccs;
  else
    if all(nccs<=1),
      nflies_per_roi(i) = 1;
    else
      nflies_per_roi(i) = 2;
    end
  end
end
roidata.nflies_per_roi = nflies_per_roi;

roifilename = fullfile(rootdatadir,experiment_name,'roiinfo.mat');
save(roifilename,'-struct','roidata');

figure(1);
clf;
axes('Position',[0,0,1,1]);
imagesc(readframe(1));
axis image;
colormap gray;
axis off;
hold on;

colors = jet(nrois)*.7;
for i = 1:nrois,
  h = drawellipse(roidata.centerx(i),roidata.centery(i),...
    0,roidata.radii(i),roidata.radii(i),'-','Color',colors(i,:));
  text(roidata.centerx(i),roidata.centery(i),sprintf('%d: %d flies',i,nflies_per_roi(i)),...
    'Color',colors(i,:),'HorizontalAlignment','center','VerticalAlignment','middle');
end

roipngname = fullfile(rootdatadir,experiment_name,'roiinfo.png');
savefig(roipngname,1,'png');


%% initialize

nframes_track = min(lastframetrack,nframes)-firstframetrack+1;
trxx = nan(2,nrois,nframes_track);
trxy = nan(2,nrois,nframes_track);
trxa = nan(2,nrois,nframes_track);
trxb = nan(2,nrois,nframes_track);
trxtheta = nan(2,nrois,nframes_track);
trxarea = nan(2,nrois,nframes_track);
mixpred = gmm(2,2,'full');
mixpred = repmat(mixpred,[1,nrois]);
predx = nan(2,nrois);
predy = nan(2,nrois);
predtheta = nan(2,nrois);
predarea = nan(2,nrois);
istouching = nan(nrois,nframes_track);


hell = nan(2,nrois);
htrx = nan(2,nrois);

gmmfit_nbadpriors = zeros(1,nrois);

%%

for t = firstframetrack:min(lastframetrack,nframes),
  
  %% backsub
  
  % read in frame
  if DEBUG >= 2,
    fprintf('Frame %d\n',t);
    tic;
  end
  im = readframe(t);
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
    
    if nflies_per_roi(roii) == 0,
      continue;
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
    obsa = nan(1,2);
    obsb = nan(1,2);
    obstheta = nan(1,2);
    obsx = nan(1,2);
    obsy = nan(1,2);
    obsarea = nan(1,2);
    obspriors = nan(1,2);
    
    cc.Area = cellfun(@numel,cc.PixelIdxList);
    ncc_bigenough = nnz(cc.Area >= minccarea);

    if DEBUG >= 2,
      fprintf('Connected components analysis took %f seconds\n',toc);
      tic;
    end

    %% only one fly
    
    if nflies_per_roi(roii) == 1,
      
      % store in position 1
      i = 1;
      if ncc_bigenough == 1,
        % use this connected component
        j = find(cc.Area >= minccarea);
        idx = cc.PixelIdxList{j};
      else
        idx = cat(1,cc.PixelIdxList{:});
      end
      [y,x] = ind2sub(size(isforebb),idx);
      w = dbkgdbb(idx);
      [mu,S] = weighted_mean_cov([x(:),y(:)],w);
      [obsa(i),obsb(i),obstheta(i)] = cov2ell(S);
      obsx(i) = mu(1);
      obsy(i) = mu(2);
      obsarea(i) = numel(cc.PixelIdxList{j});
      obspriors(i) = 1;
      
      if DEBUG >= 2,
        fprintf('Not clustering, fitting Gaussians took %f seconds\n',toc);
        tic;
      end

      istouching(roii,t-firstframetrack+1) = 0;
      order = [1,2];
      
      % store
      trxx(:,roii,iframe) = obsx;
      trxy(:,roii,iframe) = obsy;
      trxa(:,roii,iframe) = obsa;
      trxb(:,roii,iframe) = obsb;
      trxtheta(:,roii,iframe) = obstheta;
      trxarea(:,roii,iframe) = obsarea;
          
    %% fit gaussians
    
    else
    
    if ncc_bigenough == 2,
      ccidx = find(cc.Area >= minccarea);
      for i = 1:2,
        j = ccidx(i);
        [y,x] = ind2sub(size(isforebb),cc.PixelIdxList{j});
        w = dbkgdbb(cc.PixelIdxList{j});
        [mu,S] = weighted_mean_cov([x(:),y(:)],w);
        [obsa(i),obsb(i),obstheta(i)] = cov2ell(S);
        obsx(i) = mu(1);
        obsy(i) = mu(2);
        obsarea(i) = numel(cc.PixelIdxList{j});
        obspriors(i) = sum(w);
      end
      obspriors = obspriors / sum(obspriors);
      if DEBUG >= 2,
        fprintf('Not clustering, fitting Gaussians took %f seconds\n',toc);
        tic;
      end
      istouching(roii,t-firstframetrack+1) = 0;

    elseif cc.NumObjects == 0,
      fprintf('No flies detected in ROI %d, frame %d\n',roii,t);
      
    else
      %keyboard;
      % GMM clustering of all foreground pixels
      idx = cat(1,cc.PixelIdxList{:});
      [y,x] = ind2sub(size(isforebb),idx);
      w = dbkgdbb(idx);
      if t == firstframetrack,
        [mu,S,obspriors,post,nll,mixprev] = mygmm([x,y],2,...
          'Replicates',gmmem_params.nrestarts_firstframe,...
          'precision',gmmem_params.precision,...
          'MaxIters',gmmem_params.maxiters,...
          'weights',w);
      else
        [mu,S,obspriors,post,nll,mixprev] = mygmm([x,y],2,...
          'Start',mixpred(roii),...
          'precision',gmmem_params.precision,...
          'MaxIters',gmmem_params.maxiters,...
          'weights',w);
        
        % check that all went well
        if any(obspriors <= gmmem_params.min_obsprior),
          fprintf('Bad prior found, trying to reinitialize\n');
          [mu1,S1,obspriors1,post1,nll1,mixprev1] = mygmm([x,y],2,...
            'Replicates',gmmem_params.nrestarts_firstframe,...
            'precision',gmmem_params.precision,...
            'MaxIters',gmmem_params.maxiters,...
            'weights',w);
          gmmfit_nbadpriors(roii) = gmmfit_nbadpriors(roii)+1;
          if nll1 <= nll,
            fprintf('Using results from reinitialization, which improve nll by %f\n',nll-nll1);
            mu = mu1;
            S = S1;
            obspriors = obspriors1;
            post = post1;
            nll = nll1;
            mixprev = mixprev1;
          else
            fprintf('Reinitialization does not improve nll.\n');
          end
        end
        
      end
      
      for i = 1:2,
        [obsa(i),obsb(i),obstheta(i)] = cov2ell(S(:,:,i));
        obsx(i) = mu(i,1);
        obsy(i) = mu(i,2);
        obsarea(i) = sum(post(:,i));
      end
      if DEBUG >= 2,
        fprintf('GMM fitting took %f seconds\n',toc);
        tic;
      end
      istouching(roii,t-firstframetrack+1) = 1;
      
    end
   
    %% match
    
    iframe = t-firstframetrack+1;
    
    order = 1:2;
    if t > firstframetrack,
      
      besterr = inf;
      for i = 1:2,
        if i == 1,
          ordercurr = [1,2];
        else
          ordercurr = [2,1];
        end
        
        dpos2 = (predx(:,roii)'-obsx(ordercurr)).^2 + (predy(:,roii)'-obsy(ordercurr)).^2;
        dtheta = abs(modrange(predtheta(:,roii)'-obstheta(ordercurr),-pi/2,pi/2));
        darea = abs(predarea(:,roii)'-obsarea(ordercurr));
        
        errcurr = sqrt(sum(dpos2))*err_params.weightpos + ...
          sqrt(sum(dtheta.^2))*err_params.weighttheta + ...
          sqrt(sum(darea.^2))*err_params.weightarea;
        
        if errcurr < besterr,
          order = ordercurr;
          besterr = errcurr;
        end
      end
      
    end
    
    if DEBUG >= 2,
      fprintf('Matching took %f seconds\n',toc);
      tic;
    end

    end
    
    
    %% store
    trxx(:,roii,iframe) = obsx(order);
    trxy(:,roii,iframe) = obsy(order);
    trxa(:,roii,iframe) = obsa(order);
    trxb(:,roii,iframe) = obsb(order);
    trxtheta(:,roii,iframe) = obstheta(order);
    trxarea(:,roii,iframe) = obsarea(order);
        
    %% predicted position for next frame
    predarea(:,roii) = trxarea(:,roii,iframe);
    if t == firstframetrack,
      predx(:,roii) = trxx(:,roii,iframe);
      predy(:,roii) = trxy(:,roii,iframe);
      predtheta(:,roii) = trxtheta(:,roii,iframe);
    else
      predx(:,roii) = ((2-err_params.dampen_pos)*trxx(:,roii,iframe) - (1-err_params.dampen_pos)*trxx(:,roii,iframe-1));
      predy(:,roii) = ((2-err_params.dampen_pos)*trxy(:,roii,iframe) - (1-err_params.dampen_pos)*trxy(:,roii,iframe-1));
      dtheta = modrange(trxtheta(:,roii,iframe)-trxtheta(:,roii,iframe-1),-pi/2,pi/2);
      predtheta(:,roii) = (trxtheta(:,roii,iframe)+(1-err_params.dampen_theta)*dtheta);
    end
   
    % switch around priors, move toward .5, .5
    mixpred(roii).priors = (1-err_params.dampen_priors)*obspriors(order) + err_params.dampen_priors*.5;
    % set centres, covars to predicted positions
    mixpred(roii).centres = [predx(:,roii),predy(:,roii)];
    mixpred(roii).covars = axes2cov(trxa(:,roii,iframe),trxb(:,roii,iframe),predtheta(:,roii));
    
    if DEBUG >= 2,
      fprintf('Prediction took %f seconds\n',toc);
      tic;
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
  drawnow;
  
end

%% correct for bounding box of rois
for roii = 1:nrois,
  roibb = roibbs(roii,:);
  trxx(:,roii,:) = trxx(:,roii,:) + roibb(1) - 1;
  trxy(:,roii,:) = trxy(:,roii,:) + roibb(3) - 1;
end

%% resolve head/tail ambiguity

for roii = 1:nrois,

  for i = 1:2,
    x = reshape(trxx(i,roii,:),[1,nframes_track]);
    y = reshape(trxy(i,roii,:),[1,nframes_track]);
    theta = reshape(trxtheta(i,roii,:),[1,nframes_track]);
    trxtheta(i,roii,:) = choose_orientations(x,y,theta,choose_orientations_params.velocity_angle_weight,choose_orientations_params.max_velocity_angle_weight);
  end
  
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

%% assign identities based on size

assignids_nflips = nan(1,nrois);
muafit = nan(2,nrois);
sigmaafit = nan(2,nrois);
mubfit = nan(2,nrois);
sigmabfit = nan(2,nrois);
muareafit = nan(2,nrois);
sigmaareafit = nan(2,nrois);
niters_assignids_em = nan(1,nrois);
cost_assignids = nan(1,nrois);
sigmamotionfit = nan(1,nrois);
idsfit = nan(2,nrois,nframes_track);

for roii = 1:nrois,
  
  if nflies_per_roi(roii) < 2,
    continue;
  end
  
  x = reshape(trxx(:,roii,:),[2,nframes_track]);
  y = reshape(trxy(:,roii,:),[2,nframes_track]);
  a = reshape(trxa(:,roii,:),[2,nframes_track]);
  b = reshape(trxb(:,roii,:),[2,nframes_track]);
  theta = reshape(trxtheta(:,roii,:),[2,nframes_track]);
  area = a.*b.*pi;
  appearanceweight = double(~istouching(roii,:));
  iddata = cat(3,area,a,b);
  %iddata = cat(3,a,b);
  [idsfit_curr,mudatafit_curr,sigmadatafit_curr,sigmamotionfit_curr,cost_curr,niters_curr] = ...
    AssignIdentities(x,y,iddata,'vel_dampen',err_params.dampen_pos,'appearanceweight',appearanceweight);
  assignids_nflips(roii) = nnz(idsfit_curr(1,1:end-1)~=idsfit_curr(1,2:end));
  fprintf('Roi %d, flipped ids %d times\n',roii,assignids_nflips(roii));
  idx = sub2ind([2,nframes_track],idsfit_curr,repmat(1:nframes_track,[2,1]));
  trxx(:,roii,:) = x(idx);
  trxy(:,roii,:) = y(idx);
  trxa(:,roii,:) = a(idx);
  trxb(:,roii,:) = b(idx);
  trxtheta(:,roii,:) = theta(idx);
  muafit(:,roii) = mudatafit_curr(:,2);
  sigmaafit(:,roii) = sigmadatafit_curr(:,2);
  mubfit(:,roii) = mudatafit_curr(:,3);
  sigmabfit(:,roii) = sigmadatafit_curr(:,3);
  muareafit(:,roii) = mudatafit_curr(:,1);
  sigmaareafit(:,roii) = sigmadatafit_curr(:,1);
  niters_assignids_em(roii) = niters_curr;
  cost_assignids(roii) = cost_curr;
  sigmamotionfit(roii) = sigmamotionfit_curr;  
  idsfit(:,roii,:) = idsfit_curr;
end

idx = find(nflies_per_roi == 2);
areas_male = muareafit(1,idx);
areas_female = muareafit(2,idx);
area_thresh = (mean(areas_male)+mean(areas_female))/2;

idx = find(nflies_per_roi == 1);
sexperroi = cell(1,nrois);
for ii = 1:numel(idx),
  i = idx(ii);
  areacurr = trxa(1,i,:).*trxb(1,i,:)*pi;
  meanareacurr = nanmean(areacurr(:));
  fprintf('%d: %f\n',i,meanareacurr);
  if meanareacurr <= area_thresh, 
    sexperroi{i} = 'M';
  else
    sexperroi{i} = 'F';
  end
end


%% plot updated results

flycolors = {'b','r'};

for iframe = firstframetrack:nframes_track,
  
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
  end
  drawnow;
  %input('');
  
end

%% save results

timestamps = ufmf_read_timestamps(headerinfo,1,nframes);
timestamp_tracked = now;
 
% trx = struct('x',mat2cell(reshape(trxx,[2*nrois,nframes_track]),ones(1,2*nrois),nframes_track),...
%   'y',mat2cell(reshape(trxy,[2*nrois,nframes_track]),ones(1,2*nrois),nframes_track),...
%   'a',mat2cell(reshape(trxa/2,[2*nrois,nframes_track]),ones(1,2*nrois),nframes_track),...
%   'b',mat2cell(reshape(trxb/2,[2*nrois,nframes_track]),ones(1,2*nrois),nframes_track),...
%   'theta',mat2cell(reshape(trxtheta,[2*nrois,nframes_track]),ones(1,2*nrois),nframes_track));
trx = struct;
j = 1;
for i = 1:nrois,
  for jj = 1:nflies_per_roi(i),
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
    trx(j).arena = struct;
    trx(j).arena.arena_radius_mm = roidata.radii(i);
    trx(i).arena.arena_center_mm_x = roidata.centerx(i);
    trx(i).arena.arena_center_mm_y = roidata.centery(i);
    %trx(j).roipts = rois{i};
    %trx(j).roibb = roibbs(i,:);
    trx(j).moviefile = moviefile;
    trx(j).dt = diff(timestamps(trx(j).firstframe:trx(j).endframe));
    trx(j).timestamp = timestamps(trx(j).firstframe:trx(j).endframe);
    j = j + 1;
  end
  if nflies_per_roi(i) == 2,
    trx(j-2).sex = repmat({'M'},[1,nframes_track]);
    trx(j-1).sex = repmat({'F'},[1,nframes_track]);
  elseif nflies_per_roi(i) == 1,
    trx(j-1).sex = repmat(sexperroi(i),[1,nframes_track]);
  end
end    
outfilename = fullfile(rootdatadir,experiment_name,'trx.mat');
save(outfilename,'trx','timestamps','timestamp_tracked',...
  'moviefile','annfile','bg_nframes',...
  'firstframetrack','lastframetrack',...
  'bgmode','bgthresh','minccarea',...
  'gmmem_params','err_params',...
  'choose_orientations_params','bgmed','roidata');


%% track wings

expdir = fullfile(rootdatadir,experiment_name);
moviefilestr = 'movie.ufmf';
annfilestr = 'movie.ufmf.ann';
trxfilestr = 'trx.mat';
outtrxfilestr = 'wingtrx.mat';
perframedir = 'perframe';
isarena = false(imheight,imwidth);
[XGRID,YGRID] = meshgrid(1:imwidth,1:imheight);
for roii = 1:nrois,
  if nflies_per_roi(roii) == 0,
    continue;
  end
  isarena = isarena | ...
    ( ((XGRID - roidata.centerx(roii)).^2 + ...
       (YGRID - roidata.centery(roii)).^2) ...
      <= roidata.radii(i)^2 );
end

TrackWings(expdir,...
  'moviefilestr',moviefilestr,...
  'trxfilestr',trxfilestr,...
  'outtrxfilestr',outtrxfilestr,...
  'perframedir',perframedir,...
  'isarena',isarena,...
  'paramsfile','WingTrackingParameters_CourtshipBowls20121205.txt',...
  'debug',1,...
  'bgmodel',double(bgmed));
%'restart',true);
outfilename = fullfile(expdir,outtrxfilestr);
tmp = load(outfilename);
fns = {'timestamp_tracked',...
  'moviefile','annfile','bg_nframes',...
  'firstframetrack','lastframetrack',...
  'bgmode','bgthresh','minccarea',...
  'gmmem_params','err_params',...
  'choose_orientations_params','bgmed','rois'
  };
for i = 1:numel(fns),
  tmp.(fns{i}) = eval(fns{i});
end
save(outfilename,'-struct','tmp');
  

%% output tracking diagnostics

timestamp = now;
outfilename = fullfile(rootdatadir,experiment_name,'TrackTwoFliesDiagnostics.mat');
save(outfilename,...
  'timestamp','moviefile',...
  'assignids_nflips','muafit','sigmaafit',...
  'mubfit','sigmabfit','niters_assignids_em',...
  'cost_assignids','sigmamotionfit','istouching',...
  'gmmfit_nbadpriors');