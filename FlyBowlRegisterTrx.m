function trx = FlyBowlRegisterTrx(expdir,varargin)

fns_notperframe = {'id','moviename','annname','firstframe','arena','off',...
  'nframes','endframe','matname','fps','pxpermm'};

[analysis_protocol,settingsdir,registrationparamsfilestr,datalocparamsfilestr,dotemporalreg] = ...
  myparse(varargin,...
  'analysis_protocol','current',...
  'settingsdir','/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings',...
  'registrationparamsfilestr','registration_params.txt',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'dotemporalreg',true);

%% read in the data locations
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);
      
%% read in registration params

% name of parameters file
registrationparamsfile = fullfile(settingsdir,analysis_protocol,registrationparamsfilestr);
if ~exist(registrationparamsfile,'file'),
  error('Registration params file %s does not exist',registrationparamsfile);
end
% read
registration_params = ReadParams(registrationparamsfile);

%% detect registration marks

% name of annotation file
annfile = fullfile(expdir,dataloc_params.annfilestr);

% name of movie file
moviefile = fullfile(expdir,dataloc_params.moviefilestr);

% template filename should be relative to settings directory
if isfield(registration_params,'bowlMarkerType'),
  if ischar(registration_params.bowlMarkerType),
    if ~ismember(registration_params.bowlMarkerType,{'gradient'}),
      registration_params.bowlMarkerType = fullfile(settingsdir,analysis_protocol,registration_params.bowlMarkerType);
    end
  else
    % plate -> bowlmarkertype
    plateids = str2double(registration_params.bowlMarkerType(1:2:end-1));
    bowlmarkertypes = registration_params.bowlMarkerType(2:2:end);
    metadata = parseExpDir(expdir);
    plateid = str2double(metadata.plate);
    i = find(plateid == plateids,1);
    if isempty(i),
      error('bowlMarkerType not set for plate %d',plateid);
    end
    if ~ismember(bowlmarkertypes{i},{'gradient'}),
      registration_params.bowlMarkerType = fullfile(settingsdir,analysis_protocol,bowlmarkertypes{i});
    end
  end
end

% maxDistCornerFrac_BowlLabel might depend on bowl
if isfield(registration_params,'maxDistCornerFrac_BowlLabel') && ...
    numel(registration_params.maxDistCornerFrac_BowlLabel) > 1,
  plateids = registration_params.maxDistCornerFrac_BowlLabel(1:2:end-1);
  cornerfracs = registration_params.maxDistCornerFrac_BowlLabel(2:2:end);
  metadata = parseExpDir(expdir);
  plateid = str2double(metadata.plate);
  i = find(plateid == plateids,1);
  if isempty(i),
    error('maxDistCornerFrac_BowlLabel not set for plate %d',plateid);
  end
  registration_params.maxDistCornerFrac_BowlLabel = cornerfracs(i);
end

fnsignore = intersect(fieldnames(registration_params),...
  {'minFliesLoadedTime','maxFliesLoadedTime','extraBufferFliesLoadedTime'});
  
registration_params_cell = struct2paramscell(rmfield(registration_params,fnsignore));
% file to save image to
if isfield(dataloc_params,'registrationimagefilestr'),
  registration_params_cell(end+1:end+2) = {'imsavename',fullfile(expdir,dataloc_params.registrationimagefilestr)};
end
% detect
try
  registration_data = detectRegistrationMarks(registration_params_cell{:},'annName',annfile,'movieName',moviefile);
catch ME,
  error(['Error detecting registration marks:\n',getReport(ME)]);
end


%% apply spatial registration

% name of input trx mat file
ctraxfile = fullfile(expdir,dataloc_params.ctraxfilestr);

% name of movie
moviefile = fullfile(expdir,dataloc_params.moviefilestr);

% load trajectories
[trx,~,succeeded,timestamps] = load_tracks(ctraxfile,moviefile,'annname',annfile);
if ~succeeded,
  error('Could not load trajectories from file %s',ctraxfile);
end

if isempty(trx),
  error('No flies tracked.');
end

% apply to trajectories
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
  [xnose1,ynose1] = registration_data.registerfn(xnose0,ynose0);
  [xtail1,ytail1] = registration_data.registerfn(xtail0,ytail0);
  [xleft1,yleft1] = registration_data.registerfn(xleft0,yleft0);
  [xright1,yright1] = registration_data.registerfn(xright0,yright0);
  % compute the center as the mean of these four points
  x1 = (xnose1+xtail1+xleft1+xright1)/4;
  y1 = (ynose1+ytail1+yleft1+yright1)/4;
  % compute the major axis from the nose to tail distance
  a1 = sqrt( (xnose1-xtail1).^2 + (ynose1-ytail1).^2 ) / 4;
  % compute the minor axis length as the left to right distance
  b1 = sqrt( (xleft1-xright1).^2 + (yleft1-yright1).^2 ) / 4;
  % compute the orientation from the nose and tail points only
  theta1 = atan2(ynose1-ytail1,xnose1-xtail1);
  % store the registerd positions
  trx(fly).x_mm = x1;
  trx(fly).y_mm = y1;
  trx(fly).a_mm = a1;
  trx(fly).b_mm = b1;
  trx(fly).theta_mm = theta1;
  
  % add dt
  % TODO: fix timestamps after fix errors!
  if isfield(trx,'timestamps'),
    trx(fly).dt = diff(trx(fly).timestamps);
  else
    trx(fly).dt = repmat(1/trx(fly).fps,[1,trx(fly).nframes-1]);
  end

  trx(fly).fps = 1/mean(trx(fly).dt);
  trx(fly).pxpermm = 1 / registration_data.scale;
  
end

%% crop start and end of trajectories

if dotemporalreg,
  
  % how long did we actually record for?
  timestamps = timestamps - timestamps(1);
  
  % how long did we want to record for? read from config file
  configfile = dir(fullfile(expdir,dataloc_params.configfilepattern));
  if isempty(configfile),
    error('Could not find config file for experiment %s',expdir);
  end
  configfile = fullfile(expdir,configfile(1).name);
  config_params = ReadParams(configfile);
  
  % read flies loaded time -- could use metadata tree loader, but this is
  % pretty simple
  metadatafile = fullfile(expdir,dataloc_params.metadatafilestr);
  if ~exist(metadatafile,'file'),
    error('Could not find metadata file %s',metadatafile);
  end
  fid = fopen(metadatafile,'r');
  if fid < 0,
    error('Could not open metadata file %s for reading',metadatafile);
  end
  while true,
    s = fgetl(fid);
    if ~ischar(s),
      fclose(fid);
      error('Could not find "seconds_fliesloaded" in metadata file %s',metadatafile);
    end
    m = regexp(s,'seconds_fliesloaded\s*=\s*"([^"]*)"','tokens','once');
    if ~isempty(m),
      seconds_fliesloaded = str2double(m{1});
      break;
    end
  end
  fclose(fid);
  
  % how much time should we crop from the beginning?
  timeCropStart = registration_params.maxFliesLoadedTime - seconds_fliesloaded;
  if timeCropStart < 0,
    warning('Load time = %f seconds, greater than max allowed load time = %f seconds.',...
      seconds_fliesloaded,registration_params.maxFliesLoadedTime);
    timeCropStart = 0;
  end
  
  % find closest timestamp to timeCropStart
  i0 = find(timestamps >= timeCropStart,1);
  if isempty(i0),
    error('No timestamps occur after timeCropStart = %f. Cannot crop start.',timeCropStart);
  end
  if i0 > 1 && (timeCropStart-timestamps(i0-1)) < (timestamps(i0)-timeCropStart),
    i0 = i0 - 1;
  end
  
  % how long is the video currently?
  recordLengthCurr = timestamps(end)-timestamps(i0);
  
  % how long should the video be?
  recordLengthIdeal = config_params.RecordTime - registration_params.extraBufferFliesLoadedTime - ...
    (registration_params.maxFliesLoadedTime - registration_params.minFliesLoadedTime);
  
  % how much time should we crop from the end?
  if recordLengthCurr < recordLengthIdeal,
    warning('Cropped video is %f seconds long, shorter than ideal length %f seconds.',recordLengthCurr,recordLengthIdeal);
    i1 = numel(timestamps);
  else
    i1 = find(timestamps - timestamps(i0) >= recordLengthIdeal,1);
    if isempty(i1),
      warning('No timestamps occur after timeCropEnd = %f. Cannot crop end.',timestamps(i0)+recordLengthIdeal);
      i1 = numel(timestamps);
    else
      if i1 > 1 && ...
          (recordLengthIdeal - (timestamps(i1-1)-timestamps(i0))) < ...
          ((timestamps(i1)-timestamps(i0)) - recordLengthIdeal),
        i1 = i1 - 1;
      end
    end
  end
  registration_data.seconds_crop_start = timestamps(i0);
  registration_data.start_frame = i0;
  registration_data.seconds_crop_end = timestamps(end)-timestamps(i1);
  registration_data.end_frame = i1;
  
  fns = setdiff(fieldnames(trx),fns_notperframe);
  isperframe = true(1,numel(fns));
  nperfn = nan(1,numel(fns));
  if ~isempty(trx),
    
    for j = 1:numel(fns),
      fn = fns{j};
      for i = 1:numel(trx),
        if ~isnumeric(trx(i).(fn)),
          isperframe(j) = false;
          break;
        end
        if i == 1,
          ncurr = trx(i).nframes - numel(trx(i).(fn));
        else
          if trx(i).nframes - numel(trx(i).(fn)) ~= ncurr,
            isperframe(j) = false;
            break;
          end
        end
      end
      if all([trx.nframes]) == trx(1).nframes && ...
          trx(1).nframes > 1 && numel(trx(1).(fn)) == 1,
        isperframe(j) = false;
      end
      if isperframe(j),
        nperfn(j) = ncurr;
      end
    end
    
    nperfn = nperfn(isperframe);
    fns = fns(isperframe);
    ncropright = ceil(nperfn/2);
    ncropleft = nperfn - ncropright;
    
    trxdelete = false(1,numel(trx));
    for i = 1:numel(trx),
      if trx(i).firstframe > i1,
        trxdelete(i) = true;
        continue;
      end
      
      if trx(i).endframe < i0,
        trxdelete(i) = true;
        continue;
      end
      
      trx(i).nframes = min(i1,trx(i).endframe)-max(i0,trx(i).firstframe)+1;
      
      if trx(i).firstframe < i0,
        off = i0 - trx(i).firstframe;
        for j = 1:numel(fns),
          fn = fns{j};
          trx(i).(fn) = trx(i).(fn)(off+1+ncropleft(j):end);
        end
        trx(i).firstframe = i0;
      end
      
      if trx(i).endframe > i1,
        for j = 1:numel(fns),
          fn = fns{j};
          trx(i).(fn) = trx(i).(fn)(1:trx(i).nframes-nperfn(j));
        end
        trx(i).endframe = i1;
      end
      
      
      trx(i).off = -trx(i).firstframe + 1;
      
    end
    trx(trxdelete) = [];
    
  end
  
end

%% save registered trx to file

% name of output trx mat file
trxfile = fullfile(expdir,dataloc_params.trxfilestr);

try
  save(trxfile,'trx','timestamps');
catch ME
  warning('Could not save registered trx: %s',getReport(ME));
end

%% save params to mat file

registrationmatfile = fullfile(expdir,dataloc_params.registrationmatfilestr);
tmp = rmfield(registration_data,'registerfn'); %#ok<NASGU>
try
  save(registrationmatfile,'-struct','tmp');
catch ME
  warning('Could not save registered data to mat file: %s',getReport(ME));
end

%% save params to text file
registrationtxtfile = fullfile(expdir,dataloc_params.registrationtxtfilestr);
try
fid = fopen(registrationtxtfile,'w');
fnssave = {'offX','offY','offTheta','scale','bowlMarkerTheta','featureStrengths',...
  'circleCenterX','circleCenterY','circleRadius',...
  'seconds_crop_start','seconds_crop_end','start_frame','end_frame'};
fnssave = intersect(fnssave,fieldnames(registration_data));
for i = 1:numel(fnssave),
  fn = fnssave{i};
  fprintf(fid,'%s,%f\n',fn,registration_data.(fn));
end
fclose(fid);
catch ME
  warning('Could not save registered data to trx file: %s',getReport(ME));
end