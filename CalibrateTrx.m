function trx = CalibrateTrx(inTrxName,varargin)

if ~exist(inTrxName,'file'),
  error('matfile %s containing trx does not exist',inTrxName);
end

[annName,movieName,incalibrationfile,...
  outCalibrationFile,detectCalibrationParams,...
  outTrxName] = ...
  myparse(varargin,...
  'annname','',...
  'moviename','',...
  'incalibrationfile','',...
  'outcalibrationfile','',...
  'detectcalibrationparams',{},...
  'outtrxname','');

if ~isempty(incalibrationfile),
  calibrationData = load(incalibrationfile);
  calibrationData = detectCalibrationMarks('calibrationData',calibrationData);
else
   
  if ~exist(annName,'file'),
    error('annfile %s does not exist',annName);
  end
  if ~exist(movieName,'file'),
    error('moviefile %s does not exist',movieName);
  end

  % get video size
  [~,~,fid,headerinfo] = get_readframe_fcn(movieName);
  nr = headerinfo.nr;
  nc = headerinfo.nc;
  fclose(fid);
  
  % detect the calibration marks
  calibrationData = detectCalibrationMarks(detectCalibrationParams{:},...
    'annName',annName,'saveName',calibrationMarkName,...
    'nr',nr,'nc',nc);
end

% apply to trajectories
trx = load_tracks(inTrxName,movieName,'annname',annName);
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
  [xnose1,ynose1] = calibrationData.registerfn(xnose0,ynose0);
  [xtail1,ytail1] = calibrationData.registerfn(xtail0,ytail0);
  [xleft1,yleft1] = calibrationData.registerfn(xleft0,yleft0);
  [xright1,yright1] = calibrationData.registerfn(xright0,yright0);
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

if ~isempty(outCalibrationFile),
  tmp = rmfield(calibrationData,'registerfn');
  save(outCalibrationFile,'-struct',tmp);
end

if ~isempty(outTrxName),
  save(outTrxName,'trx');
end