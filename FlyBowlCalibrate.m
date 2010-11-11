function trx = FlyBowlCalibrate(expdir,protocol,varargin)

requiredParams = {...
  'calibrationMarkFileStr',...
  'movieFileStr',...
  'annFileStr',...
  'inTrxFileStr',...
  'outTrxFileStr'...
  };

% load parameters
params = FlyBowlProtocol2Params('FlyBowlCalibrate',protocol,varargin{:});
paramsRead = ismember(requiredParams,params(1,:));
if ~all(paramsRead),
  error(['The following required parameters were not read: ',sprintf('\n%s',requiredParams{~paramRead})]);
end
% add to workspace
for i = 1:size(params,2),
  eval(sprintf('%s = params{2,%d};',params{1,i},i));
end
  
% load parameters for mark detection
detectCalibrationMarksParams = ...
  FlyBowlProtocol2Params('detectCalibrationMarks',protocol,varargin{:});

% movie to calibrate
movieName = fullfile(expdir,movieFileStr);
annName = fullfile(expdir,annFileStr);
inTrxName = fullfile(expdir,inTrxFileStr);

% get video size
[readframe,nframes,fid,headerinfo] = get_readframe_fcn(movieName);
nr = headerinfo.nr;
nc = headerinfo.nc;
fclose(fid);
clear readframe;

% file to save results to
calibrationMarkName = fullfile(expdir,calibrationMarkFileStr);
outTrxName = fullfile(expdir,outTrxFileStr);

% detect the calibration marks
registerFn = detectCalibrationMarks(detectCalibrationMarksParams{:},...
  'annName',annName,'saveName',calibrationMarkName,...
  'nr',nr,'nc',nc);

% apply to trajectories
trx = load_ctrax(inTrxName,movieName,annName);
for fly = 1:length(trx),
  
  % apply transformation to 4 extremal points on the ellipse
  xnose0 = trx(fly).x + 2*trx(fly).a.*cos(trx(fly).theta);
  ynose0 = trx(fly).y + 2*trx(fly).a.*sin(trx(fly).theta);
  xtail0 = trx(fly).x - 2*trx(fly).a.*cos(trx(fly).theta);
  ytail0 = trx(fly).y - 2*trx(fly).a.*sin(trx(fly).theta);
  xleft0 = trx(fly).x + 2*trx(fly).b.*cos(trx(fly).theta-pi/2);
  yleft0 = trx(fly).y + 2*trx(fly).b.*sin(trx(fly).theta-pi/2);
  xright0 = trx(fly).x + 2*trx(fly).b.*cos(trx(fly).theta+pi/2);
  yright0 = trx(fly).y + 2*trx(fly).b.*sin(trx(fly).theta+pi/2);
  [xnose1,ynose1] = registerFn(xnose0,ynose0);
  [xtail1,ytail1] = registerFn(xtail0,ytail0);
  [xleft1,yleft1] = registerFn(xleft0,yleft0);
  [xright1,yright1] = registerFn(xright0,yright0);
  % compute the center as the mean of these four points
  x1 = (xnose1+xtail1+xleft1+xright1)/4;
  y1 = (ynose1+ytail1+yleft1+yright1)/4;
  % compute the major axis from the nose to tail distance
  a1 = sqrt( (xnose1-xtail1).^2 + (ynose1-ytail1).^2 ) / 4;
  % compute the minor axis length as the left to right distance
  b1 = sqrt( (xleft1-xright1).^2 + (yleft1-yright1).^2 ) / 4;
  % compute the orientation from the nose and tail points only
  theta1 = atan2(ynose1-ytail1,xnose1-xtail1);
  % store the calibrated positions
  trx(fly).x_mm = x1;
  trx(fly).y_mm = y1;
  trx(fly).a_mm = a1;
  trx(fly).b_mm = b1;
  trx(fly).theta_mm = theta1;
end

save(outTrxName,'trx');