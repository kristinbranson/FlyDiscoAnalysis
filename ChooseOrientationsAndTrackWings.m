function [trx,perframedata,info,wingtracking_units,trackdata] = ...
  ChooseOrientationsAndTrackWings(trx,moviefile,varargin)

min_jump_speed = 5.7256;
min_speed = 0.2264;
nframes_speed = 2;
min_wing_area = 13.0844;
max_wing_area = 22.7860;
max_ecc_confident = 0.6102;
min_ecc_factor = 0.0100;
lambda_phi = 0;
lambda_wingarea = 1;
lambda_theta = 10;

[annfile,bgmodel,isarena,bgnframes...
  wingtracking_params,...
  paramsfile,...
  min_jump_speed,...
  min_speed,...
  nframes_speed,...
  min_wing_area,...
  max_wing_area,...
  max_ecc_confident,...
  min_ecc_factor,...
  lambda_phi,...
  lambda_wingarea,...
  lambda_theta,...
  firstframe,endframe,...
  DEBUG,trackdata,...
  savefile,dofliporientation] = myparse(varargin,...
  'annfile','','bgmodel',[],'isarena',[],'bgnframes',100,...
  'wingtracking_params',[],...
  'paramsfile','',...
  'min_jump_speed',min_jump_speed,...
  'min_speed',min_speed,...
  'nframes_speed',nframes_speed,...
  'min_wing_area',min_wing_area,...
  'max_wing_area',max_wing_area,...
  'max_ecc_confident',max_ecc_confident,...
  'min_ecc_factor',min_ecc_factor,...
  'lambda_phi',lambda_phi,...
  'lambda_wingarea',lambda_wingarea,...
  'lambda_theta',lambda_theta,...
  'firstframe',1,'endframe',inf,...
  'debug',0,'trackdata',[],...
  'savefile','',...
  'dofliporientation',true);

% if trx file
if ischar(trx),
  trxfile = trx;
  trx = load_tracks(trxfile);
end

% didflip = cell(1,numel(trx));
% for fly = 1:numel(trx),
%   didflip{fly} = rand(1,trx(fly).nframes)>=.5;
%   trx(fly).theta(didflip{fly}) = modrange(trx(fly).theta(didflip{fly})+pi,-pi,pi);
% end

% if parameters in a file
if ~isempty(paramsfile),
  load(paramsfile);
end

if isempty(wingtracking_params)
  wingtracking_params = DefaultWingTrackingParams();
elseif ischar(wingtracking_params),
  [~,~,ext] = fileparts(wingtracking_params);
  switch ext,
    case '.txt'
      wingtracking_params = ReadParams(wingtracking_params);
    case '.mat'
      wingtracking_params = load(wingtracking_params);
  end
end

endframe = min(max([trx.endframe]),endframe);
  
nflies = numel(trx);

istouchinginfo = false;
iswingareainfo = false;
perframedata = [];
info = [];
wingtracking_units = [];
if lambda_wingarea > 0,
  
  if isempty(trackdata) || ~all(isfield(trackdata,{'wing_areal','wing_arear','istouching'})),
  
    if isempty(bgmodel),
      [readframe,nframes,fid] = get_readframe_fcn(moviefile);
      im = readframe(1);
      [nr,nc,~] = size(im);
      
      if ~isempty(annfile),
        [bgmodel] = read_ann(annfile,'background_center');
        bgmodel = reshape(bgmodel,[nc,nr])';
      else
        
        iscolor = size(im,3) == 3;
        ts = unique(round(linspace(1,nframes,bgnframes)));
        im = repmat(im(:,:,1),[1,1,numel(ts)]);
        for i = 1:numel(ts),
          imcurr = readframe(ts(i));
          if iscolor,
            imcurr = rgb2gray(imcolor);
          end
          im(:,:,i) = imcurr;
        end
        if ~isfloat(im),
          im = single(im);
        end
        bgmodel = median(im,3);
        
      end
      if fid > 1,
        fclose(fid);
      end
    end
    if isempty(isarena),
      [readframe,nframes,fid] = get_readframe_fcn(moviefile); %#ok<ASGLU>
      [nr,nc,~] = size(readframe(1));
      if fid > 1,
        fclose(fid);
      end
      if ~isempty(annfile),
        [isarena] = read_ann(annfile,'isarena');
        isarena = reshape(isarena,[nc,nr])' > 0;
      else
        isarena = true(nr,nc);
      end
    end
    
    % we will set to 0 wing tracking info when flies are in the same
    % connected component, so we don't care that the other flies may be in the
    % wrong orientation
    [trx0,perframedata0,info,wingtracking_units] = TrackWingsHelper(trx,moviefile,bgmodel,...
      isarena,wingtracking_params,...
      'firstframe',firstframe,...
      'endframe',endframe,...
      'debug',DEBUG);
    
    trx1 = trx;
    for fly = 1:nflies,
      trx1(fly).theta = modrange(trx(fly).theta+pi,-pi,pi);
    end
    
    [trx1,perframedata1] = TrackWingsHelper(trx1,moviefile,bgmodel,...
      isarena,wingtracking_params,...
      'firstframe',firstframe,...
      'endframe',endframe,...
      'debug',DEBUG);
    
    trackdata.wing_areal = cell(1,nflies);
    trackdata.wing_arear = cell(1,nflies);
    for fly = 1:nflies,
      trackdata.wing_areal{fly} = [perframedata1.wing_areal{fly};perframedata0.wing_areal{fly}];
      trackdata.wing_arear{fly} = [perframedata1.wing_arear{fly};perframedata0.wing_arear{fly}];
      trackdata.istouching{fly} = perframedata0.istouching{fly} | perframedata1.istouching{fly};
    end
    
    trackdata.trx0 = trx0;
    trackdata.trx1 = trx1;
    trackdata.perframedata0 = perframedata0;
    trackdata.perframedata1 = perframedata1;
    trackdata.wingtrackinginfo = info;
    trackdata.wingtracking_units = wingtracking_units;
    
  else
    
    trx0 = trackdata.trx0;
    trx1 = trackdata.trx1;
    perframedata0 = trackdata.perframedata0;
    perframedata1 = trackdata.perframedata1;
    info = trackdata.wingtrackinginfo;
    wingtracking_units = trackdata.wingtracking_units;
    
  end
  
  istouchinginfo = true;
  iswingareainfo = true;

  perframedata = perframedata0;
  trx = trx0;
  perframefns = fieldnames(perframedata0);
  
  trxwingfns = {
    'wing_anglel'
    'wing_angler'
    'xwingl'
    'ywingl'
    'xwingr'
    'ywingr'
    };
  
else
  
  trackdata = [];
  
end

if dofliporientation,
for fly = 1:nflies,
  
  t0 = max(firstframe,trx(fly).firstframe);
  t1 = min(endframe,trx(fly).endframe);
  if t1<t0,
    continue;
  end
  i0 = t0-trx(fly).firstframe + 1;
  i1 = t1-trx(fly).firstframe + 1;
  x = trx(fly).x(i0:i1);
  y = trx(fly).y(i0:i1);
  theta = trx(fly).theta(i0:i1);
  N = numel(x);
  
  % appearance cost is made up of velocity-based cost and wing-based cost
  appearancecost = zeros(2,N);
  
  % use speed to determine how confident we are in the velocity-based
  % measurement
  fil = zeros(1,nframes_speed);
  fil(1) = -1;
  fil(end) = 1;
  dx = imfilter(x,fil,'replicate');
  dy = imfilter(y,fil,'replicate');
  dt = imfilter(1:numel(x),fil,'same','corr','replicate');
  % v(i) is the speed from frame i-floor(nframes_speed/2) to i+ceil(nframes_speed/2)
  v = sqrt(dx.^2 + dy.^2)./dt;
  phi = atan2(dy,dx);
  
  dphitheta0 = abs(modrange(phi-theta,-pi,pi));
  dphitheta1 = abs(modrange(phi-theta-pi,-pi,pi));
  
  if istouchinginfo,  
    istouching = trackdata.istouching{fly}(i0:i1);
  else
    istouching = false(1,N);
  end
  
  isjumping = v>=min_jump_speed;
  
  isbad = isjumping | istouching;
  fil1 = ones(1,nframes_speed);
  isbad = imfilter(isbad,fil1,'replicate');
  
  weight_phi = max(0,min(1,(v-min_speed)/(min_jump_speed-min_speed)));
  weight_phi(isbad) = 0;
  appearancecost(1,:) = appearancecost(1,:) + dphitheta0.*weight_phi*lambda_phi;
  appearancecost(2,:) = appearancecost(2,:) + dphitheta1.*weight_phi*lambda_phi;
  
  if iswingareainfo,
    wing_area = trackdata.wing_areal{fly}(:,i0:i1) + trackdata.wing_arear{fly}(:,i0:i1);
    wing_area(:,istouching) = 0;
    
    wing_area_weight = min(1,max(0,sqrt(wing_area)-min_wing_area)/(max_wing_area-min_wing_area));
    appearancecost = appearancecost + lambda_wingarea*wing_area_weight;
  end
  
  % don't rely on change in orientation as much when the fly is circular
  % ecc_factor is >= 1 for ecc <= max_ecc_confident
  % ecc_factor is min_ecc_factor for ecc = 1
  % and interpolates linearly between (max_ecc_confident,1) and
  % (1,min_ecc_factor)
  a = trx(fly).a;
  b = trx(fly).b;
  ecc = b ./ a;
  parama = (min_ecc_factor-1)/(1-max_ecc_confident);
  paramb = min_ecc_factor - parama;
  ecc_factor = parama.*ecc + paramb;
  weight_theta = max(min_ecc_factor,min(1,ecc_factor)) .* lambda_theta;

  [newtheta,s] = choose_orientations_generic(theta,weight_theta,appearancecost); 
  trx(fly).theta(i0:i1) = newtheta;
  if isfield(trx,'theta_mm')
    absdtheta = abs(modrange(newtheta-theta,-pi,pi));
    ANGLECHANGED_THRESH = 0.1; % anything between ~eps and ~pi prob fine
    assert( isequal(absdtheta>ANGLECHANGED_THRESH,s==2) );
    tfflipped = (s==2);

    thmm = trx(fly).theta_mm(i0:i1);
    thmm(tfflipped) = modrange(thmm(tfflipped)-pi,-pi,pi);
    trx(fly).theta_mm(i0:i1) = thmm;
  end
  
  if iswingareainfo,
    idx2 = find(s==2)+i0-1;
    for i = 1:numel(perframefns),
      fn = perframefns{i};
      perframedata.(fn){fly}(idx2) = perframedata1.(fn){fly}(idx2);
    end
    for i = 1:numel(trxwingfns),
      fn = trxwingfns{i};
      trx(fly).(fn)(idx2) = trx1(fly).(fn)(idx2);
    end
    
  end  
    
end

end

if ~isempty(savefile),
  save(savefile,'trx','perframedata','info','wingtracking_units','trackdata');
end