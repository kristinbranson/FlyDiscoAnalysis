function registration = detectRegistrationMarks(varargin)

%% parse inputs

[saveName,...
  bkgdImage,movieName,annName,bkgdNSampleFrames,...
  method,...
  circleImageType,circleRLim,circleXLim,circleYLim,...
  circleImageThresh,circleCannyThresh,circleCannySigma,...
  circleNXTry,circleNYTry,circleNRTry,...
  crossFilterRadius,nRotations,minTemplateFeatureStrength,...
  diskFilterRadius,minFeatureStrengthLow,minFeatureStrengthHigh,...
  minDistCenterFrac,maxDistCenterFrac,...
  maxDistCornerFrac_BowlLabel,...
  nRegistrationPoints,featureRadius,...
  maxDThetaMate,...
  pairDist_mm,...
  circleRadius_mm,...
  bowlMarkerPairTheta_true,...
  maxDThetaBowlMarkerPair,...
  markerPairAngle_true,...
  nBowlMarkers,...
  bowlMarkerType,...
  DEBUG,...
  registration,...
  nr,nc,...
  imsavename,...
  hfig,figpos,...
  useNormXCorr,...
  doTransposeImage,...
  ledindicator,...
  regXY] = ...
  myparse(varargin,...
  'saveName','',...
  'bkgdImage',[],...
  'movieName','',...
  'annName','',...
  'bkgdNSampleFrames',10,...
  'method','normcorr',...
  'circleImageType','canny',...
  'circleRLim',[.475,.51],...
  'circleXLim',[.475,.525],...
  'circleYLim',[.475,.525],...
  'circleImageThresh',1,...
  'circleCannyThresh',[],...
  'circleCannySigma',[],...
  'circleNXTry',50,...
  'circleNYTry',50,...
  'circleNRTry',50,...
  'crossFilterRadius',9,...
  'nRotations',20,...
  'minTemplateFeatureStrength',.92,...
  'diskFilterRadius',11,...
  'minFeatureStrengthLow',20,...
  'minFeatureStrengthHigh',30,...
  'minDistCenterFrac',.5,...
  'maxDistCenterFrac',.57,...
  'maxDistCornerFrac_BowlLabel',.17,...
  'nRegistrationPoints',8,...
  'featureRadius',25,...
  'maxDThetaMate',10*pi/180,...
  'pairDist_mm',133,...
  'circleRadius_mm',63.5,...
  'bowlMarkerPairTheta_true',-3*pi/4,...
  'maxDThetaBowlMarkerPair',pi/12,...
  'markerPairAngle_true',pi/6,...
  'nBowlMarkers',1,...
  'bowlMarkerType','gradient',...
  'debug',false,...
  'registrationData',[],...
  'nr',[],'nc',[],...
  'imsavename','',...
  'hfig',1,'figpos',[10,10,1600,800],...
  'useNormXCorr',false,...
  'doTransposeImage',false,...
  'ledindicator',false,...
  'regXY',[]);
%annName = [] ;

iscircle = ismember(method,{'circle','circle_manual'});

%% return register function

if ~isempty(registration),
  registration.registerfn = @(x,y) register(x,y,registration.offX,...
    registration.offY,registration.offTheta,registration.scale);
  registration.affine = affineTransform(registration.offX,...
    registration.offY,registration.offTheta,registration.scale);
  return;
end

%% get background image

isBkgdImage = ~isempty(bkgdImage);

if ~isBkgdImage && ~isempty(annName),
  % try reading 
  [bkgdImage,bkgdMed,bkgdMean,bg_algorithm,movie_height,movie_width] = ...
    read_ann(annName,'background_center',...
    'background_median','background_mean','bg_algorithm',...
    'movie_height','movie_width');
  if isempty(bkgdImage),
    if strcmpi(bg_algorithm,'median'),
      bkgdImage = bkgdMed;
    else
      bkgdImage = bkgdMean;
    end
  end
  if isempty(bkgdImage),
    error('Could not read background center from ann file');
  end
  if ~isempty(movie_height),
    nr = movie_height;
  end
  if ~isempty(movie_width),
    nc = movie_width;
  end
  if isempty(nr) || isempty(nc),
    if exist(movieName,'file'),
      [readframe,~,fid,~] = get_readframe_fcn(movieName);
      im = readframe(1);
      nr = size(im,1);
      nc = size(im,2);
      fclose(fid);
    else
      error('Shape of movie could not be read from ann file of movie file');
    end
  end
  if doTransposeImage,
    bkgdImage = reshape(bkgdImage,[nc,nr])';
  else
    bkgdImage = reshape(bkgdImage,[nr,nc]);
  end
  isBkgdImage = true;
end
if ~isBkgdImage,
  if isempty(movieName),
    error('Either bkgdImage, annName, or movieName must be input. Ann file must contain background center parameter');
  end
  if ~exist(movieName,'file'),
    error('Movie file %s does not exist',movieName);
  end
  % open the movie
  [readframe,~,fid,headerinfo] = get_readframe_fcn(movieName);
  
  % take the median
  if isfield(headerinfo, 'nmeans') && headerinfo.nmeans > 1,
    meanims = ufmf_read_mean(headerinfo,'meani',2:headerinfo.nmeans);
  else
    sampleframes = unique(round(linspace(1,headerinfo.nframes,bkgdNSampleFrames)));
    meanims = repmat(double(readframe(1)),[1,1,1,bkgdNSampleFrames]);
    for i = 2:bkgdNSampleFrames,
      meanims(:,:,:,i) = double(readframe(sampleframes(i)));
    end
  end
  meanims = double(meanims);
  bkgdImage = median(meanims,4);
  if fid>0 ,
      fclose(fid);
  end
end

[nr,nc,ncolors] = size(bkgdImage); %#ok<NASGU>
r = min(nr,nc);


%% feature detection filtering

% compute gradient magnitude image
gradI = [diff(bkgdImage,1,1).^2;zeros(1,nc)] + [diff(bkgdImage,1,2).^2,zeros(nr,1)];
% filter with uniform filter of input radius
fil = fspecial('disk',diskFilterRadius);
gradfilI = imfilter(gradI,fil,0,'same');

if strcmpi(method,'normcorr'),
  
  % make a cross filter
  fil0 = zeros(crossFilterRadius*2+1);
  fil0(crossFilterRadius+1,:) = 1;
  fil0(:,crossFilterRadius+1) = 1;
  
  % steer
  fils = cell(1,nRotations);
  thetas = linspace(0,90,nRotations+1);
  thetas = thetas(1:end-1);
  for i = 1:nRotations,
    fils{i} = 1-2*imrotate(fil0,thetas(i),'bilinear','loose');
  end

  % compute normalized maximum correlation
  filI1 = -inf(nr,nc);
  for i = 1:nRotations,
    filI1 = max(filI1,imfilter(bkgdImage,fils{i},'replicate') ./ ...
                      imfilter(bkgdImage,ones(size(fils{i})),'replicate'));
  end
  
elseif strcmpi(method,'gradient'),
  
  filI1 = gradfilI;
  
elseif strcmpi(method,'circle'),
  
  % hough circle transform to detect 
  if strcmpi(circleImageType,'raw_whiteedge'),
    circleim = bkgdImage >= circleImageThresh;
  elseif strcmpi(circleImageType,'raw_blackedge'),
    circleim = bkgdImage <= circleImageThresh;
  elseif strcmpi(circleImageType,'grad'),
    circleim = sqrt(gradI) >= circleImageThresh;
  elseif strcmpi(circleImageType,'canny'),
    circleim = bkgdImage;
  else
    error('Unknown circleImageType %s',circleImageType);
  end
  binedgesa = linspace(circleXLim(1)*nc,circleXLim(2)*nc,circleNXTry+1);
  bincentersb = linspace(circleYLim(1)*nr,circleYLim(2)*nr,circleNYTry);
  bincentersr = linspace(circleRLim(1)*min(nc,nr),circleRLim(2)*min(nc,nr),circleNRTry);
  [circleRadius,circleCenterX,circleCenterY,featureStrengths,circleDetectParams] = ...
    detectcircles(circleim,...
    'cannythresh',circleCannyThresh,'cannysigma',circleCannySigma,...
    'binedgesa',binedgesa,'bincentersb',bincentersb,'bincentersr',bincentersr,...
    'maxncircles',1,'doedgedetect',strcmpi(circleImageType,'canny'));

elseif strcmpi(method,'circle_manual'),
  
  circleim = bkgdImage;
  hfig = figure;
  imagesc(bkgdImage,[0,255]); axis image;
  [circleCenterX,circleCenterY,circleRadius] = fitcircle_manual(hfig);
  if ishandle(hfig),
    delete(hfig);
  end
  featureStrengths = nan;
  
else
  
  error('Unknown method for registration mark detection: %s',method);
  
end


%% find bowl label

% compute distance to corners
if nBowlMarkers > 0,
  [xGrid,yGrid] = meshgrid(1:nc,1:nr);
  [dxGrid,dyGrid] = meshgrid(-featureRadius:featureRadius,-featureRadius:featureRadius);
  distCorner = inf(nr,nc);
  corners = [1,1,nc,nc;1,nr,nr,1];
  for i = 1:size(corners,2),
    distCorner = min(distCorner, sqrt( (xGrid-corners(1,i)).^2 + (yGrid-corners(2,i)).^2 ));
  end
  
  % threshold max distance to some corner
  switch bowlMarkerType,
    case 'gradient',
      filI4 = gradfilI;
      methodcurr = 'grad2';
    otherwise,
      bowlMarkerTemplate = im2double(imread(bowlMarkerType));
      h = fspecial('gaussian',[5,5],1);
      bowlMarkerTemplate =  imfilter(bowlMarkerTemplate,h);
      bkgdImage = imfilter(bkgdImage,h);
      
      if useNormXCorr,
        bowlMarkerTemplate = 255*bowlMarkerTemplate;
        sz1 = size(bowlMarkerTemplate);
        offr0 = floor(sz1(1)/2);
        offr1 = sz1(1)-offr0;
        offc0 = floor(sz1(2)/2);
        offc1 = sz1(2)-offc0;
      else %is what runs in ledma
        bowlMarkerTemplate = bowlMarkerTemplate - min(bowlMarkerTemplate(:));
        bowlMarkerTemplate = bowlMarkerTemplate / max(bowlMarkerTemplate(:));
        bowlMarkerTemplate = 2*bowlMarkerTemplate-1;
      end
      bowlfils = cell(1,2*nRotations);
      thetas = linspace(0,180,2*nRotations+1);
      thetas = thetas(1:end-1);
      for i = 1:2*nRotations,
        bowlfils{i} = imrotate(bowlMarkerTemplate,thetas(i),'bilinear','loose');
      end
      % compute normalized maximum correlation
      filI4 = -inf(nr,nc);
      for i = 1:2*nRotations,
        if useNormXCorr,
          tmp = normxcorr2(bowlMarkerTemplate,bkgdImage);
          tmp = tmp(offr0+1:end-offr1+1,offc0+1:end-offc1+1);
          tmp(1:sz1(1),:) = 0; tmp(end-sz1(1)+1:end,:) = 0;
          tmp(:,1:sz1(2)) = 0; tmp(:,end-sz1(2)+1:end) = 0;
          filI4 = max(filI4,tmp);
        else
          filI4 = max(filI4,imfilter(bkgdImage,bowlfils{i},'replicate') ./ ...
            imfilter(bkgdImage,ones(size(bowlfils{i})),'replicate'));
        end
      end
      methodcurr = 'template';
  end

  bowlMarkerIm = filI4; %#ok<NASGU>

  maxDistCorner_BowlLabel = maxDistCornerFrac_BowlLabel * r;
  filI4(distCorner > maxDistCorner_BowlLabel) = -inf;
  % make corner with registration mark -inf
  if ledindicator
    filI4 = neginfOutDetection(regXY(2),regXY(1),filI4,maxDistCorner_BowlLabel);
%     filI4 = zeroOutDetection(regXY(2),regXY(1),filI4);
  end
  
  bowlMarkerPoints = nan(2,nBowlMarkers);
  for i = 1:nBowlMarkers,
    % find maximum
    [success,x,y,featureStrength] = getNextFeaturePoint(filI4,methodcurr);
    if success,
      bowlMarkerPoints(:,i) = [x;y];
      filI4 = zeroOutDetection(bowlMarkerPoints(1,i),bowlMarkerPoints(2,i),filI4);
    else
      error('Could not detect bowl marker %d',i);
    end
  end
end

%

%% find registration points

if ~iscircle

  % compute distance from center of image
  minDistCenter = minDistCenterFrac * r;
  maxDistCenter = maxDistCenterFrac * r;
  distCenter = sqrt((xGrid-nc/2).^2 + (yGrid-nr/2).^2);
  
  % threshold distance from center
  filI2 = filI1;
  filI2(distCenter < minDistCenter | distCenter > maxDistCenter) = 0;
  
  % zero out bowl marker
  for i = 1:nBowlMarkerPoints,
    filI2 = zeroOutDetection(bowlMarkerPoints(1,i),bowlMarkerPoints(2,i),filI2);
  end
  
  registrationPoints = [];
  featureStrengths = [];
  filI3 = filI2;
  for i = 1:nRegistrationPoints,
    
    [success,x,y,featureStrength,filI3] = getNextFeaturePoint(filI3,method);
    if ~success,
      break;
    end
    registrationPoints(:,i) = [x;y]; %#ok<AGROW>
    featureStrengths(i) = featureStrength; %#ok<AGROW>
    
  end
  
  if isempty(registrationPoints),
    error('No registration points detected');
  end
  
end

%% origin is the average of all the registration points
if iscircle,
  originX = circleCenterX;
  originY = circleCenterY;
else
  isleft = registrationPoints(1,:) <= nc/2;
  nLeftRegistrationPoints = nnz(isleft);
  nRightRegistrationPoints = nnz(~isleft);
  istop = registrationPoints(2,:) >= nr/2;
  nTopRegistrationPoints = nnz(istop);
  nBottomRegistrationPoints = nnz(~istop);
  
  originX = (sum(registrationPoints(1,isleft))/nLeftRegistrationPoints + ...
    sum(registrationPoints(1,~isleft))/nRightRegistrationPoints)/2;
  originY = (sum(registrationPoints(2,istop))/nTopRegistrationPoints + ...
    sum(registrationPoints(2,~istop))/nBottomRegistrationPoints)/2;
end

offX = -originX;
offY = -originY;

%% sort registration points counterclockwise from bowl label

if nBowlMarkers > 0 && success,
  bowlMarkerPoint = nanmean(bowlMarkerPoints,2);
  bowlMarkerTheta = atan2(bowlMarkerPoint(2)-originY,bowlMarkerPoint(1)-originX);
else
  bowlMarkerTheta = atan2(1-originY,1-originX);
end
if ~iscircle
  theta = atan2(registrationPoints(2,:)-originY,registrationPoints(1,:)-originX);
  % offset from bowlMarkerTheta
  dtheta = mod((theta - bowlMarkerTheta),2*pi);
  [dtheta,order] = sort(dtheta);
  registrationPoints = registrationPoints(:,order);
  featureStrengths = featureStrengths(order);
end

%% find the rotation

% if nBowlMarkers > 0 && success,
% 
%   if iscircle
%     % registration marker
%     bowlMarkerPairTheta = bowlMarkerTheta;
%   else
%     % take the average of the pair of points around the bowl marker
%     d1 = mod(dtheta(1)+pi,2*pi)-pi;
%     d2 = mod(dtheta(end)+pi,2*pi)-pi;
%     inrange = d1 > 0 && d1 <= maxDThetaBowlMarkerPair && ...
%       d2 < 0 && -d2 <= maxDThetaBowlMarkerPair;
%     if inrange,
%       x = (registrationPoints(1,1)+registrationPoints(1,end))/2;
%       y = (registrationPoints(2,1)+registrationPoints(2,end))/2;
%       bowlMarkerPairTheta = atan2(y-originY,x-originX);
%     else
%       bowlMarkerPairTheta = bowlMarkerTheta;
%     end
%   end
%   
%   offTheta = mod(bowlMarkerPairTheta_true-bowlMarkerPairTheta+pi,2*pi)-pi;
% else
%   offTheta = 0;
% end

offTheta = 0;
%% find scale

if iscircle,
  scale = circleRadius_mm / circleRadius;
else
  nPairs = 0;
  d = 0;
  for i = 1:size(registrationPoints,2),
    thetaCurr = dtheta(i);
    if thetaCurr > pi,
      break;
    end
    [dThetaMate,j] = min(abs(thetaCurr+pi-dtheta));
    isMate = dThetaMate <= maxDThetaMate;
    if isMate,
      d = d + sqrt(diff(registrationPoints(1,[i,j])).^2+diff(registrationPoints(2,[i,j])).^2);
      nPairs = nPairs + 1;
    end
  end
  meanDistPair_px = d / nPairs;
  scale = pairDist_mm / meanDistPair_px;
end

registerfn = @(x,y) register(x,y,offX,offY,offTheta,scale);

affine = affineTransform(offX,offY,offTheta,scale);

registration = struct('offX',offX,...
  'offY',offY,...
  'offTheta',offTheta,...
  'scale',scale,...
  'bowlMarkerTheta',bowlMarkerTheta,...
  'bkgdImage',bkgdImage,...
  'featureStrengths',featureStrengths,...
  'affine',affine);

registration.bowlMarkerPoints = bowlMarkerPoints;

if iscircle,
  registration.circleCenterX = circleCenterX;
  registration.circleCenterY = circleCenterY;
  registration.circleRadius = circleRadius;
  if strcmpi(method,'circle'),
    registration.circleDetectParams = circleDetectParams;
  end
else
  registration.registrationPoints = registrationPoints;
  registration.bowlMarkerPoints = bowlMarkerPoints;
  registration.bowlMarkerPoint = bowlMarkerPoint;
  registration.dtheta = dtheta;
  registration.nPairs = nPairs;
  registration.nLeftRegistrationPoints = nLeftRegistrationPoints;
  registration.nRightRegistrationPoints = nRightRegistrationPoints;
  registration.nTopRegistrationPoints = nTopRegistrationPoints;
  registration.nBottomRegistrationPoints = nBottomRegistrationPoints;
end
if ~isempty(saveName),
  save(saveName,'-struct','registration');
end
registration.registerfn = registerfn;

%% create images illustrating fitting

if ~isempty(imsavename) && ~ledindicator,
  
  figure(hfig);
  clf(hfig);
  set(hfig,'Units','pixels','Position',figpos);
  nimsplot = 2;
  hax = createsubplots(1,nimsplot,.025,hfig);
  
  % plot background image in jet colormap
  imagesc(bkgdImage,'parent',hax(1),[0,255]);
  axis(hax(1),'image');
  axis(hax(1),'xy');
  hold(hax(1),'on');
  
  % plot detected registration points if not circle mode
  if ~iscircle,
    plot(hax(1),registrationPoints(1,:),registrationPoints(2,:),'ks');
    plot(hax(1),registrationPoints(1,:),registrationPoints(2,:),'wd');
  end

  % plot detected/labeled circle if is circle mode
  if iscircle,
    tmp = linspace(0,2*pi,ceil(2*pi*circleRadius));
    plot(hax(1),circleCenterX + circleRadius*cos(tmp),...
      circleCenterY + circleRadius*sin(tmp),...
      'w-','linewidth',2);
    plot(hax(1),circleCenterX + circleRadius*cos(tmp),...
      circleCenterY + circleRadius*sin(tmp),...
      'k--','linewidth',2);
    plot(hax(1),circleCenterX,circleCenterY,'ks','markerfacecolor','k');
    plot(hax(1),circleCenterX,circleCenterY,'wd');
  end
  
  % plot bowl marker(s)
  if nBowlMarkers > 0,
    plot(hax(1),bowlMarkerPoints(1,:),bowlMarkerPoints(2,:),'ks','markersize',12);
    plot(hax(1),bowlMarkerPoints(1,:),bowlMarkerPoints(2,:),'kx','markersize',12);
  end
  
  % plot axes
  xangle = 0;
  yangle = pi/2;
  if iscircle,
    l = circleRadius_mm/2 / scale;
  else
    l = pairDist_mm/4 / scale;
  end
  quiver(hax(1),originX,originY,cos(xangle-offTheta)*l,sin(xangle-offTheta)*l,0,'k-');
  quiver(hax(1),originX,originY,cos(xangle-offTheta)*l,sin(xangle-offTheta)*l,0,'w--');
  text(originX+cos(xangle-offTheta)*l,originY+sin(xangle-offTheta)*l,'x','parent',hax(1));
  quiver(hax(1),originX,originY,cos(yangle-offTheta)*l,sin(yangle-offTheta)*l,0,'w-');
  quiver(hax(1),originX,originY,cos(yangle-offTheta)*l,sin(yangle-offTheta)*l,0,'k--');
  text(originX+cos(yangle-offTheta)*l,originY+sin(yangle-offTheta)*l,'y','parent',hax(1));
  
  title(hax(1),'Input background image');

  % transform image
  A = affineTransform(offX,offY,offTheta,scale);
  T = maketform('affine',A);
  % transformations of corners
  %udata = [1;nc]; vdata = [1;nr];
  %outbounds = findbounds(T,[udata,vdata]);
  %xdata = outbounds(:,1);
  %ydata = outbounds(:,2);
  [imreg,xdata,ydata] = imtransform(bkgdImage,T,'XYScale',scale,'FillValues',nan);
  
  % draw registered image
  imagesc(xdata,ydata,imreg,'parent',hax(2),[0,255]);
  hold(hax(2),'on');
    
  % draw angle to bowl marker
  if iscircle,
    l = circleRadius_mm;
  else
    l = pairDist_mm;
  end  
  plot(hax(2),[0,cos(bowlMarkerPairTheta_true)]*l,...
    [0,sin(bowlMarkerPairTheta_true)]*l,'k-o','markerfacecolor','k');
  plot(hax(2),[0,cos(bowlMarkerPairTheta_true)]*l,...
    [0,sin(bowlMarkerPairTheta_true)]*l,'w:');
  text(cos(bowlMarkerPairTheta_true)*l/2,...
    sin(bowlMarkerPairTheta_true)*l/2,sprintf('%.2fmm',l));
  text(cos(bowlMarkerPairTheta_true)*l/5,...
    sin(bowlMarkerPairTheta_true)*l/5,sprintf('%.1fdeg',bowlMarkerPairTheta_true*180/pi));
  
  % draw axes
  if iscircle,
    l = circleRadius_mm/2;
  else
    l = pairDist_mm/4;
  end
  quiver(hax(2),0,0,cos(xangle)*l,sin(xangle)*l,0,'k-');
  quiver(hax(2),0,0,cos(xangle)*l,sin(xangle)*l,0,'w--');
  text(cos(xangle)*l,sin(xangle)*l,'x','parent',hax(2));
  quiver(hax(2),0,0,cos(yangle)*l,sin(yangle)*l,0,'w-');
  quiver(hax(2),0,0,cos(yangle)*l,sin(yangle)*l,0,'k--');
  text(cos(yangle)*l,sin(yangle)*l,'y','parent',hax(2));

  axis(hax(2),'image');
  axis(hax(2),'xy');
  title(hax(2),'Registered bkgd image');
  colormap(hfig,'jet');

  try
    if exist(imsavename,'file'),
      try
        delete(imsavename);
      catch ME,
        warning('Could not delete file %s:\n %s',imsavename,getReport(ME));
      end
    end
    set(hfig,'Units','pixels','Position',figpos);
    save2png(imsavename,hfig);
  catch ME,
    warning('Could not save to %s:\n%s',imsavename,getReport(ME));
  end
  
  if isdeployed && ishandle(hfig),
    close(hfig);
  end
  
end
%% create image illustrating LED detection
if ~isempty(imsavename) && ledindicator,
  figpos = [10,10,800,800];
  figure(hfig);
  clf(hfig);
  set(hfig,'Units','pixels','Position',figpos);
  nimsplot = 1;
  hax = createsubplots(1,nimsplot,.025,hfig);
   
  % plot background image in jet colormap
  imagesc(bkgdImage,'parent',hax(1),[0,1]);
  axis(hax(1),'image');
  axis(hax(1),'xy');
  hold(hax(1),'on');

  % plot detected registration points if not circle mode
  if ~iscircle,
    plot(hax(1),registrationPoints(1,:),registrationPoints(2,:),'ks');
    plot(hax(1),registrationPoints(1,:),registrationPoints(2,:),'wd');
  end

  % plot detected/labeled circle if is circle mode
  if iscircle,
    tmp = linspace(0,2*pi,ceil(2*pi*circleRadius));
    plot(hax(1),circleCenterX + circleRadius*cos(tmp),...
      circleCenterY + circleRadius*sin(tmp),...
      'w-','linewidth',2);
    plot(hax(1),circleCenterX + circleRadius*cos(tmp),...
      circleCenterY + circleRadius*sin(tmp),...
      'k--','linewidth',2);
    plot(hax(1),circleCenterX,circleCenterY,'ks','markerfacecolor','k');
    plot(hax(1),circleCenterX,circleCenterY,'wd');
  end
  
  % plot bowl marker(s)
  if nBowlMarkers > 0,
    plot(hax(1),bowlMarkerPoints(1,:),bowlMarkerPoints(2,:),'ks','markersize',12);
    plot(hax(1),bowlMarkerPoints(1,:),bowlMarkerPoints(2,:),'kx','markersize',12);
  end

 title(hax(1),'Maxvalue image and LED detection');
 colormap(hfig,'jet');

  try
    if exist(imsavename,'file'),
      try
        delete(imsavename);
      catch ME,
        warning('Could not delete file %s:\n %s',imsavename,getReport(ME));
      end
    end
    set(hfig,'Units','pixels','Position',figpos);
    save2png(imsavename,hfig);
  catch ME,
    warning('Could not save to %s:\n%s',imsavename,getReport(ME));
  end
  
  if isdeployed && ishandle(hfig),
    close(hfig);
  end
  
end
 
 
%%

if DEBUG,
  hfig = figure;
  set(hfig,'Position',[20,20,1500,600]);
  clf;
  if iscircle,
    nsubplots = 3;
  else
    nsubplots = 4;
  end
  hax = createsubplots(1,nsubplots,.05);
  %hax(1) = subplot(1,nsubplots,1);
  axes(hax(1)); %#ok<MAXES>
  imagesc(bkgdImage);
  hold on;
  if iscircle,
    l = circleRadius_mm/2 / scale;
  else
    l = pairDist_mm/4 / scale;
  end
  xangle = 0;
  yangle = pi/2;
  quiver(originX,originY,cos(xangle-offTheta)*l,sin(xangle-offTheta)*l,0,'k--');
  text(originX+cos(xangle-offTheta)*l,originY+sin(xangle-offTheta)*l,'x');
  quiver(originX,originY,cos(yangle-offTheta)*l,sin(yangle-offTheta)*l,0,'k-');
  text(originX+cos(yangle-offTheta)*l,originY+sin(yangle-offTheta)*l,'y');
  if iscircle,
    tmp = linspace(0,2*pi,50);
    plot(circleCenterX + circleRadius*cos(tmp),circleCenterY + circleRadius*sin(tmp),'k-');
  else
    plot(registrationPoints(1,:),registrationPoints(2,:),'ks');
  end
  if nBowlMarkers > 0,
    plot(bowlMarkerPoints(1,:),bowlMarkerPoints(2,:),'mo');
  end
  axis image xy;
  axes(hax(2)); %#ok<MAXES>
  %hax(2) = subplot(1,nsubplots,2);
  if iscircle,
    imagesc(bkgdImage);
  else
    image(cat(3,min(1,bkgdImage/255 + double(filI1 >= minFeatureStrengthLow)*.3), repmat(bkgdImage/255,[1,1,2])));
  end
  hold on;
  quiver(originX,originY,cos(xangle-offTheta)*l,sin(xangle-offTheta)*l,0,'k--');
  text(originX+cos(xangle-offTheta)*l,originY+sin(xangle-offTheta)*l,'x');
  quiver(originX,originY,cos(yangle-offTheta)*l,sin(yangle-offTheta)*l,0,'k-');
  text(originX+cos(yangle-offTheta)*l,originY+sin(yangle-offTheta)*l,'y');
  if iscircle,
    tmp = linspace(0,2*pi,50);
    plot(circleCenterX + circleRadius*cos(tmp),circleCenterY + circleRadius*sin(tmp),'k-');
  else
    scatter(registrationPoints(1,:),registrationPoints(2,:),50,1:size(registrationPoints,2),'s');
  end
  if nBowlMarkers > 0,
    plot(bowlMarkerPoints(1,:),bowlMarkerPoints(2,:),'mo');
  end
  axis image xy;
  axes(hax(3)); %#ok<MAXES>
  %hax(3) = subplot(1,nsubplots,3);
  if iscircle,
    imagesc(circleim);
  else
    imfilI1 = colormap_image(filI1);
    image(imfilI1);
  end
  hold on;
  quiver(originX,originY,cos(xangle-offTheta)*l,sin(xangle-offTheta)*l,0,'w--');
  text(originX+cos(xangle-offTheta)*l,originY+sin(xangle-offTheta)*l,'x','color','w');
  quiver(originX,originY,cos(yangle-offTheta)*l,sin(yangle-offTheta)*l,0,'w-');
  text(originX+cos(yangle-offTheta)*l,originY+sin(yangle-offTheta)*l,'y','color','w');
  if iscircle,
    tmp = linspace(0,2*pi,50);
    plot(circleCenterX + circleRadius*cos(tmp),circleCenterY + circleRadius*sin(tmp),'k-');
  else
    scatter(registrationPoints(1,:),registrationPoints(2,:),50,1:size(registrationPoints,2),'s');
  end
  if nBowlMarkers > 0,
    plot(bowlMarkerPoints(1,:),bowlMarkerPoints(2,:),'mo');
  end
  axis image xy;
  linkaxes(hax(1:3));
  if ~iscircle
    axes(hax(4)); %#ok<MAXES>
    %subplot(1,nsubplots,4);
    tmp = linspace(0,2*pi,100);
    plot(cos(tmp)*pairDist_mm/2,sin(tmp)*pairDist_mm/2,'b');
    hold on;
    tmp = pi/2:pi/2:2*pi;
    plot([cos(bowlMarkerPairTheta_true+markerPairAngle_true/2+tmp)*pairDist_mm/2;zeros(size(tmp))],...
      [sin(bowlMarkerPairTheta_true+markerPairAngle_true/2+tmp)*pairDist_mm/2;zeros(size(tmp))],'g-');
    plot([cos(bowlMarkerPairTheta_true-markerPairAngle_true/2+tmp)*pairDist_mm/2;zeros(size(tmp))],...
      [sin(bowlMarkerPairTheta_true-markerPairAngle_true/2+tmp)*pairDist_mm/2;zeros(size(tmp))],'g-');
    plot([0,cos(bowlMarkerPairTheta_true)*pairDist_mm/2],[0,sin(bowlMarkerPairTheta_true)*pairDist_mm/2],'m-');
    [x,y] = registerfn(registrationPoints(1,:),registrationPoints(2,:));
    plot(x,y,'ks');
    [x,y] = registerfn(bowlMarkerPoint(1),bowlMarkerPoint(2));
    plot(x,y,'mo');
    quiver(0,0,cos(xangle)*pairDist_mm/4,sin(xangle)*pairDist_mm/4,0,'k');
    text(cos(xangle)*pairDist_mm/4,sin(xangle)*pairDist_mm/4,'x');
    quiver(0,0,cos(yangle)*pairDist_mm/4,sin(yangle)*pairDist_mm/4,0,'k');
    text(cos(yangle)*pairDist_mm/4,sin(yangle)*pairDist_mm/4,'y');
    quiver(0,0,0,pairDist_mm/4,0,'k');
    axis equal;
    axisalmosttight;
  end
  
end

%%

function [success,x,y,featureStrength,filI] = getNextFeaturePoint(filI,methodcurr)

  success = false;

  % find the next strongest feature
  [featureStrength,j] = max(filI(:));
  [y,x] = ind2sub([nr,nc],j);
  
  % make sure it meets threshold
  if strcmpi(methodcurr,'grad2'),
    if featureStrength < minFeatureStrengthHigh,
      return;
    end
  else
    if featureStrength < minTemplateFeatureStrength,
      return;
    end
  end

  % subpixel accuracy: 
  
  if strcmpi(methodcurr,'grad2'),
    % take box around point and compute weighted average of feature strength
    [box] = padgrab(filI,0,y-featureRadius,y+featureRadius,x-featureRadius,x+featureRadius);
    box = double(box > minFeatureStrengthLow);
    Z = sum(box(:));
    dx = sum(box(:).*dxGrid(:))/Z;
    dy = sum(box(:).*dyGrid(:))/Z;
    x = x + dx;
    y = y + dy;
  end
  
  filI = zeroOutDetection(x,y,filI);
  success = true;

end

  function filI = zeroOutDetection(x,y,filI)
    
    % zero out region around feature
    i1 = max(1,round(x)-featureRadius);
    i2 = min(nc,round(x)+featureRadius);
    j1 = max(1,round(y)-featureRadius);
    j2 = min(nr,round(y)+featureRadius);
    filI(j1:j2,i1:i2) = 0;
    
  end

 function filI = neginfOutDetection(x,y,filI,rd)
    
    % zero out region around feature
    i1 = max(1,round(round(x)-rd));  % ALT 2021-03-05: Added outer round() to eliminate warning, should preserve behavior
    i2 = min(nc,round(round(x)+rd));
    j1 = max(1,round(round(y)-rd));
    j2 = min(nr,round(round(y)+rd));
    filI(j1:j2,i1:i2) = -inf;  
    
  end

  function [x,y] = register(x,y,offX,offY,offTheta,scale)
    sz = size(x);
    if numel(sz) ~= numel(size(y)) || ~all(sz == size(y)),
      error('Size of x and y must match');
    end
    costheta = cos(offTheta); sintheta = sin(offTheta);
    X = [x(:)'+offX;y(:)'+offY];
    X = [costheta,-sintheta;sintheta,costheta] * X * scale;
    x = reshape(X(1,:),sz);
    y = reshape(X(2,:),sz);    
  end

  function A = affineTransform(offX,offY,offTheta,scale)
    costheta = cos(offTheta); sintheta = sin(offTheta);
    A = [1 0 0; 0 1 0; offX offY 1] * ...
      [costheta sintheta 0; -sintheta costheta 0; 0 0 1] * ...
      [scale 0 0; 0 scale 0; 0 0 1];
  end
end