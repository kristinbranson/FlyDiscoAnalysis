%% set up path

%rootdatadir = '/groups/branson/bransonlab/projects/olympiad/WingTrackedData';
addpath ../misc;
addpath ../filehandling;
addpath ../../Ctrax/matlab/netlab;
%addpath /groups/branson/home/bransonk/olympiad_bransonlab/FlyBowlAnalysis;
rootdatadir = '/groups/branson/bransonlab/projects/FlyBowl/Courtship_plates_movies';
experiment_names = {'CantonS_TrpA_Rig2Plate17BowlC_20121210T150843'
  'CantonS_TrpA_Rig2Plate17BowlD_20121210T150854'
  'CantonS_TrpA_Rig2Plate17BowlA_20121210T150418'
  'CantonS_TrpA_Rig2Plate17BowlB_20121210T150427'};
settingsdir = '/groups/branson/bransonlab/projects/CourtshipBowls/CourtshipBowlAnalysis/settings';
analysis_protocol = 'current';

%% parameters

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

roifilestr = 'roiinfo.mat';

DEBUG = 0;

%% loop

for expi = 1:numel(experiment_names),
%expi = 1;
experiment_name = experiment_names{expi};

%% read stuff in

moviefile = fullfile(rootdatadir,experiment_name,'movie.ufmf');
annfile = fullfile(rootdatadir,experiment_name,'movie.ufmf.ann');

[readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefile);
imheight = headerinfo.nr;
imwidth = headerinfo.nc;

%% compute background model

if ~isempty(annfile) && exist(annfile,'file'),
  bgmed = read_ann(annfile,'background_center');
  bgmed = uint8(round(reshape(bgmed,[imwidth,imheight])'));
else
  fprintf('Computing background model...\n');
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
roidata = struct;
[roidata.centerx,roidata.centery,roidata.radii,roidata.detectscores] = ...
  detectCourtshipBowlROIs(bgmed,roiparams);



% 
% [metadata,success] = parseExpDir(experiment_name);
% if success,
%   rigbowl = [metadata.rig,metadata.bowl];
%   success = false;
%   for i = 1:expi-1,
%     [metadata0,success0] = parseExpDir(experiment_names{i});
%     if ~success0,
%       continue;
%     end
%     rigbowl0 = [metadata0.rig,metadata0.bowl];
%     if strcmp(rigbowl,rigbowl0),
%       expi_prev = i;
%       success = true;
%     end
%   end
% end
% 
% if success,
%   roifilename = fullfile(rootdatadir,experiment_names{expi_prev},'roiinfo.mat');
% else
%   roifilename = fullfile(rootdatadir,experiment_names{expi},'roiinfo.mat');
% end
% 
% while true,
%   
%   hfig = 1;
%   figure(hfig);
%   clf;
%   imagesc(readframe(round((1+nframes)/2)));
%   axis image;
%   colormap gray;
%   hax = gca;
%   
%   
%   if ~exist(roifilename,'file'),
%     
%     title('Click to add next region of interest. Double-click on the polygon to continue.');
%     
%     idxroi = zeros(size(bgmed));
%     rois = {};
%     inrois = {};
%     roibbs = zeros(0,4);
%     [XGRID,YGRID] = meshgrid(1:imwidth,1:imheight);
%     
%     while true,
%       
%       h = impoly(hax,'Closed',true);
%       position = wait(h);
%       hold on;
%       hplot = plot(position(:,1),position(:,2),'r.-');
%       res = questdlg('Choose next action:','ROI next action','Add another ROI','Finish','Undo','Add another ROI');
%       switch res,
%         case {'Finish','Add another ROI'},
%           rois{end+1} = position; %#ok<SAGROW>
%           bb = [floor(max(1,floor(min(position(:,1))))),ceil(min(imwidth,(max(position(:,1))))),...
%             floor(max(1,floor(min(position(:,2))))),ceil(min(imheight,(max(position(:,2)))))];
%           roibbs(end+1,:) = bb; %#ok<SAGROW>
%           in = inpolygon(XGRID,YGRID,position(:,1),position(:,2));
%           idxroi(in) = numel(rois);
%           inrois{end+1} = in(bb(3):bb(4),bb(1):bb(2)); %#ok<SAGROW>
%           if strcmp(res,'Finish'),
%             break;
%           end
%         case 'Undo',
%           if ishandle(hplot),
%             delete(hplot);
%           end
%       end
%     end
%     
%     nrois = numel(rois);
%     
%     save(roifilename,'rois','nrois','roibbs','idxroi','inrois');
%     
%     break;
%     
%     
%   else
%     
%     load(roifilename);
%     hold on;
%     for i = 1:nrois,
%       plot(rois{i}(:,1),rois{i}(:,2),'r.-');
%     end
%     
%     res = questdlg('Happy with these rois? or enter new ones?','Happy?','Happy','New','Happy');
%     
%     if strcmpi(res,'New'),
%       roifilename = fullfile(rootdatadir,experiment_names{expi},'roiinfo.mat');
%       continue;
%     end
%     
%     break;
%     
%   end
%   
% end


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

    
    %% fit gaussians
    
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

%% save results

timestamps = ufmf_read_timestamps(headerinfo,1,nframes);
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
    trx(j).obsarea = reshape(trxarea(jj,i,:),[1,nframes_track]);
  end
  j = sub2ind([2,nrois],1,i);
  trx(j).sex = repmat({'M'},[1,nframes_track]);
  j = sub2ind([2,nrois],2,i);
  trx(j).sex = repmat({'F'},[1,nframes_track]);
end    
outfilename = fullfile(rootdatadir,experiment_name,'tmptrx.mat');
save(outfilename,'trx','timestamps','timestamp_tracked',...
  'moviefile','annfile','bg_nframes',...
  'firstframetrack','lastframetrack',...
  'bgmode','bgthresh','minccarea',...
  'gmmem_params','err_params',...
  'choose_orientations_params','bgmed','rois');




%% track wings

wingparamsfile = 'WingTrackingParameters_CourtshipBowls20121205.txt';

expdir = fullfile(rootdatadir,experiment_name);
moviefilestr = 'movie.ufmf';
annfilestr = 'movie.ufmf.ann';
trxfilestr = 'tmptrx.mat';
outtrxfilestr = 'tmpwingtrx.mat';
perframedir = 'perframe';
isarena = false(imheight,imwidth);
[XGRID,YGRID] = meshgrid(1:imwidth,1:imheight);
for roii = 1:nrois,
  isarena = isarena | inpolygon(XGRID,YGRID,rois{roii}(:,1),rois{roii}(:,2));
end

TrackWings(expdir,...
  'moviefilestr',moviefilestr,...
  'trxfilestr',trxfilestr,...
  'outtrxfilestr',outtrxfilestr,...
  'perframedir',perframedir,...
  'isarena',isarena,...
  'paramsfile',wingparamsfile,...
  'debug',DEBUG,...
  'bgmodel',double(bgmed));
%'restart',true);

wingperframedata = struct;
perframefns = {'nwingsdetected','wing_areal','wing_arear','wing_trough_angle'};

for i = 1:numel(perframefns),
  fn = perframefns{i};
  tmp = load(fullfile(expdir,perframedir,[fn,'.mat']));
  wingperframedata.(fn) = reshape(cell2mat(tmp.data'),[2,nrois,nframes_track]);
end

wingtrxdata = struct;
tmp = load(fullfile(expdir,outtrxfilestr));
trxfns = {'wing_anglel','wing_angler','xwingl','xwingr','ywingl','ywingr'};
for i = 1:numel(trxfns),
  fn = trxfns{i};
  wingtrxdata.(fn) = reshape(cat(1,tmp.trx.(fn)),[2,nrois,nframes_track]);
end

%% assign identities based on area of wing pixels

assignids_nflips = nan(1,nrois);
muwingareafit = nan(2,nrois);
sigmawingareafit = nan(2,nrois);
niters_assignids_em = nan(1,nrois);
cost_assignids = nan(1,nrois);
sigmamotionfit = nan(1,nrois);
idsfit = nan(2,nrois,nframes_track);

for roii = 1:nrois,
  x = reshape(trxx(:,roii,:),[2,nframes_track]);
  y = reshape(trxy(:,roii,:),[2,nframes_track]);
  wingarea = reshape(wingperframedata.wing_areal(:,roii,:)+wingperframedata.wing_arear(:,roii,:),[2,nframes_track]);
  appearanceweight = double(~istouching(roii,:));
%   
%   doflipid = rand(1,nframes_track) > .5;
%   idsperm = [double(doflipid)+1;2-double(doflipid)];
%   idx = sub2ind([2,nframes_track],idsperm,repmat(1:nframes_track,[2,1]));
%   x = x(idx);
%   y = y(idx);
%   wingarea = wingarea(idx);
  
  [idsfit_curr,mudatafit_curr,sigmadatafit_curr,sigmamotionfit_curr,cost_curr,niters_curr] = ...
    AssignIdentities(x,y,wingarea,'vel_dampen',err_params.dampen_pos,'appearanceweight',appearanceweight);
  assignids_nflips(roii) = nnz(idsfit_curr(1,1:end-1)~=idsfit_curr(1,2:end));
  fprintf('Roi %d, flipped ids %d times\n',roii,assignids_nflips(roii));
  idx = sub2ind([2,nrois,nframes_track],idsfit_curr,repmat(roii,[2,nframes_track]),repmat(1:nframes_track,[2,1]));
  trxx(:,roii,:) = trxx(idx);
  trxarea(:,roii,:) = trxarea(idx);
  trxy(:,roii,:) = trxy(idx);
  trxa(:,roii,:) = trxa(idx);
  trxb(:,roii,:) = trxb(idx);
  idsfit(:,roii,:) = idsfit_curr;
  trxtheta(:,roii,:) = trxtheta(idx);
  muwingareafit(:,roii) = mudatafit_curr;
  sigmawingareafit(:,roii) = sigmadatafit_curr;
  niters_assignids_em(roii) = niters_curr;
  cost_assignids(roii) = cost_curr;
  sigmamotionfit(roii) = sigmamotionfit_curr;  
  
  for i = 1:numel(perframefns),
    fn = perframefns{i};
    wingperframedata.(fn)(:,roii,:) = wingperframedata.(fn)(idx);
  end
  for i = 1:numel(trxfns),
    fn = trxfns{i};
    wingtrxdata.(fn)(:,roii,:) = wingtrxdata.(fn)(idx);
  end
  
end

%% plot updated results

flycolors = {'b','r'};

for iframe = 1:nframes_track,    
  
  %try
  
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
  
  %catch ME
  %  warning(getReport(ME));
  %  break;
  %end
  
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
    trx(j).obsarea = reshape(trxarea(jj,i,:),[1,nframes_track]);
    trx(j).sex = repmat({'M'},[1,nframes_track]);
  end
  j = sub2ind([2,nrois],1,i);
  trx(j).wingsclipped = true(1,nframes_track);
  j = sub2ind([2,nrois],2,i);
  trx(j).wingsclipped = false(1,nframes_track);
end    
outfilename = fullfile(rootdatadir,experiment_name,'trx.mat');
save(outfilename,'trx','timestamps','timestamp_tracked',...
  'moviefile','annfile','bg_nframes',...
  'firstframetrack','lastframetrack',...
  'bgmode','bgthresh','minccarea',...
  'gmmem_params','err_params',...
  'choose_orientations_params','bgmed','rois');


%% save wing tracking results

for i = 1:nrois,
  for jj = 1:2,
    j = sub2ind([2,nrois],jj,i);
    for k = 1:numel(trxfns),
      fn = trxfns{k};
      trx(j).(fn) = reshape(wingtrxdata.(fn)(jj,i,:),[1,nframes_track]);
    end
  end
end  
outfilename = fullfile(rootdatadir,experiment_name,'wingtrx.mat');
save(outfilename,'trx','timestamps','timestamp_tracked',...
  'moviefile','annfile','bg_nframes',...
  'firstframetrack','lastframetrack',...
  'bgmode','bgthresh','minccarea',...
  'gmmem_params','err_params',...
  'choose_orientations_params','bgmed','rois');

for k = 1:numel(perframefns),
  fn = perframefns{k};
  data = cell(1,nrois*2);
  for i = 1:nrois,
    for jj = 1:2,
      j = sub2ind([2,nrois],jj,i);
      data{j} = reshape(wingperframedata.(fn)(jj,i,:),[1,nframes_track]);
    end
  end
  save('-append',fullfile(expdir,perframedir,[fn,'.mat']),'data');
end

%% output tracking diagnostics

timestamp = now;
outfilename = fullfile(rootdatadir,experiment_name,'TrackTwoFliesDiagnostics.mat');
save(outfilename,...
  'timestamp','moviefile',...
  'assignids_nflips','idsfit','muwingareafit','sigmawingareafit',...
  'niters_assignids_em',...
  'cost_assignids','sigmamotionfit','istouching',...
  'gmmfit_nbadpriors');

%% end loop

fclose(fid);

end