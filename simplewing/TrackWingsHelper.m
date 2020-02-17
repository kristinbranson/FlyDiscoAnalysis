function [trx,perframedata,info,units] = TrackWingsHelper(trx,moviefile,bgmodel,isarena,params,varargin)

info.trackwings_version = '0.1';
info.trackwings_timestamp = datestr(now,'yyyymmddTHHMMSS');

%% parse parameters

tmpfilename = sprintf('TmpResultsTrackWings_%s.mat',datestr(now,'yyyymmddTHHMMSSPFFF'));

[firstframe,DEBUG,restart,tmpfilename,framestrack,perframedata,endframe] = ...
  myparse(varargin,...
  'firstframe',1,...
  'debug',0,...
  'restart',false,...
  'tmpfilename',tmpfilename,...
  'framestrack',[],...
  'perframedata',[],...
  'endframe',inf);

fprintf('TrackWings temporary results saved to file %s\n',tmpfilename);

debugdata = struct;
debugdata.DEBUG = DEBUG;

% choose histogram bins for wing pixel angles for fitting wings
params.edges_dthetawing = linspace(-params.max_wingpx_angle,params.max_wingpx_angle,params.nbins_dthetawing+1);
params.centers_dthetawing = (params.edges_dthetawing(1:end-1)+params.edges_dthetawing(2:end))/2;
params.wing_peak_min_frac = 1/params.nbins_dthetawing*params.wing_peak_min_frac_factor;

% morphology structural elements
params.se_dilate_body = strel('disk',params.radius_dilate_body);
params.se_open_wing = strel('disk',params.radius_open_wing);

% for sub-bin accuracy in fitting the wings
params.subbin_x = (-params.wing_radius_quadfit_bins:params.wing_radius_quadfit_bins)';    

% units for everything
units = struct(...
  'nwingsdetected',parseunits('unit'),...
  'wing_areal',parseunits('px^2'),...
  'wing_arear',parseunits('px^2'),...
  'wing_trough_angle',parseunits('rad'));

%% open movie

[readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefile); %#ok<NASGU,ASGLU>
[nr,nc,~] = size(readframe(1));
nflies = numel(trx);
npx = nr*nc; %#ok<NASGU>

%% allocate

[XGRID,YGRID] = meshgrid(1:nc,1:nr);

if isempty(perframedata),
  perframedata = struct;
  perframedata.nwingsdetected = cell(1,nflies);
  perframedata.wing_areal = cell(1,nflies);
  perframedata.wing_arear = cell(1,nflies);
  perframedata.wing_trough_angle = cell(1,nflies);
  perframedata.istouching = cell(1,nflies);
  
  for fly = 1:nflies,
    trx(fly).wing_anglel = nan(1,trx(fly).nframes);
    trx(fly).wing_angler = nan(1,trx(fly).nframes);
    perframedata.nwingsdetected{fly} = nan(1,trx(fly).nframes);
    perframedata.wing_areal{fly} = nan(1,trx(fly).nframes);
    perframedata.wing_arear{fly} = nan(1,trx(fly).nframes);
    perframedata.wing_trough_angle{fly} = nan(1,trx(fly).nframes);
    perframedata.istouching{fly} = false(1,trx(fly).nframes);
  end
end

% trajectories for the current frame
trxcurr = struct(...
  'x',cell(1,nflies),...
  'y',cell(1,nflies),...
  'a',cell(1,nflies),...
  'b',cell(1,nflies),...
  'theta',cell(1,nflies),...
  'firstframe',cell(1,nflies),...
  'endframe',cell(1,nflies),...
  'nframes',cell(1,nflies),...
  'off',cell(1,nflies)...
  );
  

%% initialize debug plots

if debugdata.DEBUG,
  debugdata.colors = hsv(nflies);
  debugdata.colors = debugdata.colors(randperm(nflies),:);
  debugdata.hims = nan(1,4);
  if debugdata.DEBUG > 1,
    debugdata.hfig = 1;
    figure(debugdata.hfig);
    clf;
    debugdata.hax = createsubplots(1,3,.01);
  end
  drawnow;
end

%% start tracking

wingtrxprev = [];

if ischar(restart),
  load(restart,'trx','perframedata','t','wingtrxprev');
  fprintf('Restarting tracking at frame %d...\n',t); %#ok<NODEF>
  startframe = t;
else
  startframe = firstframe;
end

if isempty(framestrack),
  framestrack = max(startframe,min([trx.firstframe])):min(endframe,max([trx.endframe]));
end

%for t = round(linspace(max(firstframe,min([trx.firstframe])),max([trx.endframe]),50)),
%for t = max(startframe,min([trx.firstframe])):max([trx.endframe]),

nframestracked = 0;
for t = framestrack(:)',
  
  nframestracked = nframestracked + 1;

  if mod(nframestracked,30) == 0,
    fprintf('Frame %d\n',t);
    drawnow;
  end

  for fly = 1:nflies,
     if t < trx(fly).firstframe || t > trx(fly).endframe,
      trxcurr(fly).x = [];
      trxcurr(fly).y = [];
      trxcurr(fly).a = [];
      trxcurr(fly).b = [];
      trxcurr(fly).theta = [];
     else
       i = trx(fly).off+t;
       trxcurr(fly).x = trx(fly).x(i);
       trxcurr(fly).y = trx(fly).y(i);
       trxcurr(fly).a = trx(fly).a(i);
       trxcurr(fly).b = trx(fly).b(i);
       trxcurr(fly).theta = trx(fly).theta(i);
    end
  end
  
  % fit this frame
  im = double(readframe(t));
  if debugdata.DEBUG,
    debugdata.im = im;
    debugdata.t = t;
  end
  [wingtrxcurr,debugdata,istouchingcurr] = TrackWingsOneFrame(im,bgmodel,isarena,trxcurr,wingtrxprev,params,XGRID,YGRID,debugdata);
  
  % store results
  for fly = 1:nflies,
    
    if t < trx(fly).firstframe || t > trx(fly).endframe,
      continue;
    end
    i = trx(fly).off+t;
    
    trx(fly).wing_anglel(i) = wingtrxcurr(fly).wing_anglel;
    trx(fly).wing_angler(i) = wingtrxcurr(fly).wing_angler;
    
    perframedata.nwingsdetected{fly}(i) = wingtrxcurr(fly).nwingsdetected;
    perframedata.wing_areal{fly}(i) = wingtrxcurr(fly).wing_areal;
    perframedata.wing_arear{fly}(i) = wingtrxcurr(fly).wing_arear;
    perframedata.wing_trough_angle{fly}(i) = wingtrxcurr(fly).wing_trough_angle;
    perframedata.istouching{fly}(i) = istouchingcurr(fly);
    
  end
  
  if ~isdeployed && mod(t,5000) == 0,
    save(tmpfilename,'trx','perframedata','t','wingtrxprev');
  end

  wingtrxprev = wingtrxcurr;

end

fprintf('Tracking complete.\n');

%% add in x, y positions for plotting

for fly = 1:nflies,    
  trx(fly).xwingl = trx(fly).x + 4*trx(fly).a.*cos(trx(fly).theta+ pi+trx(fly).wing_anglel);
  trx(fly).ywingl = trx(fly).y + 4*trx(fly).a.*sin(trx(fly).theta+ pi+trx(fly).wing_anglel);
  trx(fly).xwingr = trx(fly).x + 4*trx(fly).a.*cos(trx(fly).theta+ pi+trx(fly).wing_angler);
  trx(fly).ywingr = trx(fly).y + 4*trx(fly).a.*sin(trx(fly).theta+ pi+trx(fly).wing_angler);
end

%% clean up

% remove temporary file
if exist(tmpfilename,'file'),
  try
    delete(tmpfilename);
  catch ME,
    warning('Could not delete tmp file: %s',getReport(ME));
  end
end
   
if isdeployed && debugdata.DEBUG > 0,
  close all;
end
fclose(fid);