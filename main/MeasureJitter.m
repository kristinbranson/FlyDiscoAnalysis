function [dx,dy,xs,ys] = MeasureJitter(expdir,varargin)

[analysis_protocol,settingsdir,datalocparamsfilestr,...
  nsamples,framespersample,markerrad,croprad,endoff] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir',default_settings_folder_path(),...
  'datalocparamsfilestr','dataloc_params.txt',...
  'nsamples',10,'nframespersample',30,...
  'markerrad',50,'croprad',75,'endoff',0);


%% location of data
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

moviefile = fullfile(expdir,dataloc_params.moviefilestr);
trxfile = fullfile(expdir,dataloc_params.trxfilestr);

%% read parameters

regdata = load(fullfile(expdir,dataloc_params.registrationmatfilestr));

%% detect registration markers

[readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviefile);
nframes = nframes - endoff;

if isfield(regdata,'start_frame'),
  T0 = regdata.start_frame;
else
  T0 = 1;
end
if isfield(regdata,'end_frame'),
  T1 = regdata.end_frame;
else
  T1 = nframes;
end
nframes = T1-T0+1;

if nframes/framespersample > nsamples,
  t0s = unique(round(linspace(T0,T1-framespersample,nsamples)));
  t1s = t0s+framespersample-1;
else
  t0s = T0:framespersample:T1-framespersample+1;
  t1s = t0s+framespersample-1;
end
nsamples = numel(t0s);

% crop out registration mark from a frame
t = round(nframes/2);
im = readframe(t);
imsz = size(im);

markx = regdata.bowlMarkerPoints(2);
marky = regdata.bowlMarkerPoints(1);
markx0 = markx-markerrad;
markx1 = markx+markerrad;
if markx1 > imsz(2),
  markx1 = imsz(2);
  markx = markx1-2*markerrad;
elseif markx < 1,
  markx = 1;
  markx1 = markx+2*markerrad;
end
marky0 = marky-markerrad;
marky1 = marky+markerrad;
if marky1 > imsz(1),
  marky1 = imsz(1);
  marky0 = marky1-2*markerrad;
elseif marky0 < 1,
  marky0 = 1;
  marky1 = marky0+2*markerrad;
end
mark = double(im(marky0:marky1,markx0:markx1));

cropx0 = markx-croprad;
cropx1 = markx+croprad;
if cropx1 > imsz(2),
  cropx1 = imsz(2);
  cropx0 = cropx1-2*croprad;
elseif cropx0 < 1,
  cropx0 = 1;
  cropx1 = cropx0+2*croprad;
end
cropy0 = marky-croprad;
cropy1 = marky+croprad;
if cropy1 > imsz(1),
  cropy1 = imsz(1);
  cropy0 = cropy1-2*croprad;
elseif cropy0 < 1,
  cropy0 = 1;
  cropy1 = cropy0+2*croprad;
end

%cropsz = [cropy1-cropy0+1,cropx1-cropx0+1];
z = ones(2*markerrad+1);

xs = nan(1,nframes);
ys = nan(1,nframes);
for samplei = 1:nsamples,
  %fprintf('Sample %d / %d\n',samplei,nsamples);
  t0 = t0s(samplei);
  t1 = t1s(samplei);
  for t = t0:t1,
    im = readframe(t);
    imcrop = double(im(cropy0:cropy1,cropx0:cropx1));
    fil = conv2(imcrop,mark,'valid') ./ conv2(imcrop,z,'valid');
    [~,idx] = max(fil(:));
    [y,x] = ind2sub(size(fil),idx);
    ys(t) = y;
    xs(t) = x;
  end
end

dx = max(xs)-min(xs);
dy = max(ys)-min(ys);

fclose(fid);

if false,
  
  %%
  analysis_protocol = '20210531_flybubble_LED_AR_20210819';
  expdir = '/groups/branson/home/bransonk/behavioranalysis/code/MABe2022/data/nochr_TrpA91B01_Unknown_RigB_20201216T152556';
  MeasureJitter(expdir,'analysis_protocol',analysis_protocol);
  
  %%
end
