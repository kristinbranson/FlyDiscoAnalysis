function FlyDiscoRegisterTrx(expdir, varargin)

%
% Parse arguments
%
[analysis_protocol,settingsdir,registrationparamsfilestr,datalocparamsfilestr,dotemporalreg,dotemporaltruncation] = ...
  myparse(varargin,...
  'analysis_protocol','current_bubble',...
  'settingsdir',default_settings_folder_path(),...
  'registrationparamsfilestr','registration_params.txt',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'dotemporalreg',false,'dotemporaltruncation',false);


%
% read in the data locations
%
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);


%
% read in registration params
%

% name of parameters file
registrationparamsfile = fullfile(settingsdir,analysis_protocol,registrationparamsfilestr);
if ~exist(registrationparamsfile,'file'),
  error('Registration params file %s does not exist',registrationparamsfile);
end

% read the file contents into a struct
registration_params = ReadParams(registrationparamsfile);

% Use registration_params to determine final dotemporalreg, dotemporaltruncation
[dotemporalreg, dotemporaltruncation] = determineTemporalStuff(registration_params, dotemporalreg, dotemporaltruncation) ;


%
% Detect registration marks
%

% Load the background file output by FlyTracker
if isfield(dataloc_params,'flytrackerbgstr')
  flytrackerbgfile = fullfile(expdir,dataloc_params.flytrackerbgstr);
  load(flytrackerbgfile,'bg');
  bg_mean = 255*bg.bg_mean;
else
  % Note: We really do want to error if this is missing.
  % We've been bitten by this not being what we thought it was.
  error('dataloc_params is missing field flytrackerbgstr, which is required');
end

% Load metadata
metadata = collect_metadata(expdir, dataloc_params.metadatafilestr) ;

% Tweak registration_params.bowlMarkerType, if present
if isfield(registration_params,'bowlMarkerType'),
  registration_params.bowlMarkerType = determineBowlMarkerType(registration_params.bowlMarkerType, metadata, settingsdir, analysis_protocol) ;
end

% Tweak registration_params.maxDistCornerFrac_BowlLabel, if present
if isfield(registration_params,'maxDistCornerFrac_BowlLabel') ,
  registration_params.maxDistCornerFrac_BowlLabel = determineMaxDistCornerFracBowlLabel(registration_params.maxDistCornerFrac_BowlLabel, metadata) ;
end

% Call the core registration mark detection routine
registration_params_cell = marshallRegistrationParamsForDetectRegistrationMarks(registration_params, expdir, dataloc_params) ;
registration_data = detectRegistrationMarks(registration_params_cell{:},'bkgdImage',bg_mean,'useNormXCorr',true);
fprintf('Detected registration marks.\n');


%
% Apply spatial registration to trajectories
%

% Get name of input trx mat file
ctraxfile = fullfile(expdir,dataloc_params.ctraxfilestr);

% Get name of movie
moviefile = fullfile(expdir,dataloc_params.moviefilestr);

% Load trajectories
[trx,~,succeeded,timestamps] = load_tracks(ctraxfile,moviefile,'annname','');
if ~succeeded,
  error('Could not load trajectories from file %s',ctraxfile);
end

% Postprocess trajectories to remove nans from flytracker outputs
nids0 = numel(trx);
if ~isfield(registration_params,'maxFlyTrackerNanInterpFrames'),
  fprintf('maxFlyTrackerNanInterpFrames not set in registration_params, using default value.\n');
  args = {};
else
  args = {'maxFlyTrackerNanInterpFrames',registration_params.maxFlyTrackerNanInterpFrames};
end
[trx,ninterpframes,newid2oldid] = PostprocessFlyTracker(trx,args{:});
fprintf('Removed nans from tracker output.\n');
fprintf('Number of nans interpolated through: %d frames\n',ninterpframes);
fprintf('Number of identities was %d, now %d\n',nids0,numel(trx));

% Store some metadata about nan-removal in registration_data
registration_data.flytracker_nnanframes = ninterpframes;
registration_data.flytracker_nids0 = nids0;

% Compute the median frame interval, if called for
if registration_params.usemediandt,
  meddt = medianIntervalFromTimestamps(timestamps) ;
else
  meddt = [] ;  % Should not be accessed if registration_params.usemediandt is false
end

% If there are zero tracks, something is wrong
if isempty(trx),
  error('No flies tracked.');
end

% Apply spatial registration to trajectories
trx = appendPhysicalUnitFieldsToTrx(trx, registration_params.usemediandt, meddt) ;
fprintf('Applied spatial registration.\n') ;

%
% Handle all the stuff that is specific to optogenetic experiments
%
registration_params = handleOptogeneticExperimentRegistration(registration_params) ;


%%% Crop start and end of trajectories based on fly loaded time
% Need this 'global' in several places below  
fns_notperframe = {'id','moviename','annname','firstframe','arena','off',...
                   'nframes','endframe','matname','fps','pxpermm'};
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
  ndelete = 0;
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
    newid2oldid(trxdelete) = [];
    ndelete = nnz(trxdelete);
  end  
  fprintf('Applied temporal registration, deleted %d trajectories.\n',ndelete);  
else
  fprintf('NOT applying temporal registration.\n');  
end

%%% truncate end of movie based on value from registration params 
if dotemporaltruncation    
  ndelete = 0;
  
  headerInfo_local = ufmf_read_header(moviefile);
  if ~isfield(headerInfo_local,'timestamps'),
    error('No field timestamps in UFMF header');
  end
  timestamps_header = headerInfo_local.timestamps;    %??? couldn't figure out where timestamps come from in load_tracks for movie_JAABA/trx.mat
    
  % how long is the video
  recordLengthCurr = timestamps_header(end);
  
  % how long should the video be?
  recordLengthIdeal = registration_params.doTemporalTruncation;
  
  % how much time should we crop from the end?
  if recordLengthCurr < recordLengthIdeal,
    warning('Cropped video is %f seconds long, shorter than ideal length %f seconds.',recordLengthCurr,recordLengthIdeal);
    i1 = numel(timestamps_header);
  else
    i1 = find(timestamps_header >= recordLengthIdeal,1);
    if isempty(i1),
      warning('No timestamps occur after timeCropEnd = %f. Cannot crop end.',timestamps(i0)+recordLengthIdeal); %??? why warn and not error? 
      i1 = numel(timestamps_header);
    else
      %??? couldn't figure out what this is checking for
      if i1 > 1 && ...
          (recordLengthIdeal - timestamps_header(i1-1)) < ...
          (timestamps_header(i1) - recordLengthIdeal), 
        i1 = i1 - 1;
      end
    end
  end

  %not cropping start
  i0 = 1;
  registration_data.seconds_crop_start = 0;
  registration_data.start_frame = i0;
  registration_data.seconds_crop_end = timestamps_header(end)-timestamps_header(i1);
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
        % check if fn has perframe data 
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
    newid2oldid(trxdelete) = [];
    ndelete = nnz(trxdelete);    
  end  
  fprintf('Applied temporal truncation. Data cropped to %f seconds. Deleted %d trajectories.\n',timestamps_header(i1),ndelete);
else
  fprintf('NOT applying temporal truncation.\n')
end

%%% save registered trx to file
% name of output trx mat file
trxfile = fullfile(expdir,dataloc_params.trxfilestr);
didsave = false;
try
  if exist(trxfile,'file'),
    delete(trxfile);
  end
  save(trxfile,'trx','timestamps'); %??? which timestamps should this be? Shouldn't it be truncated to match registration data? 
  didsave = true;
catch ME
  warning('Could not save registered trx: %s',getReport(ME));
end
if didsave,
  fprintf('Saved registered trx to file %s\n',trxfile);
else
  fprintf('Could not save registered trx:\n%s\n',getReport(ME));
end

%%% save params to mat file
registration_data.newid2oldid = newid2oldid;
registration_data.flytracker_nidsnew = numel(trx);
registrationmatfile = fullfile(expdir,dataloc_params.registrationmatfilestr);
tmp = rmfield(registration_data,'registerfn'); 
didsave = false;
try
  if exist(registrationmatfile,'file'),
    delete(registrationmatfile);
  end
  save(registrationmatfile,'-struct','tmp');
  didsave = true;
catch ME
  warning('Could not save registered data to mat file: %s',getReport(ME));
end
if didsave,
  fprintf('Saved registration data to file %s\n',registrationmatfile);
else
  fprintf('Could not save registration data to mat file:\n%s\n',getReport(ME));
end

%%% save params to text file
registrationtxtfile = fullfile(expdir,dataloc_params.registrationtxtfilestr);
didsave = false;
try
  if exist(registrationtxtfile,'file'),
    delete(registrationtxtfile);
  end
  fid = fopen(registrationtxtfile,'w');
  
  fnssave = {'offX','offY','offTheta','scale','bowlMarkerTheta','featureStrengths',...
    'circleCenterX','circleCenterY','circleRadius',...
    'seconds_crop_start','seconds_crop_end','start_frame','end_frame',...
    'flytracker_nnanframes','flytracker_nids0','flytracker_nidsnew'};
  
  fnssave = intersect(fnssave,fieldnames(registration_data));
  for i = 1:numel(fnssave),
    fn = fnssave{i};
    fprintf(fid,'%s,%f\n',fn,registration_data.(fn));
  end
  
  if isfield(registration_params,'OptogeneticExp')
    if registration_params.OptogeneticExp
      fprintf(fid,'%s,%d\n','ledX',registration_data.ledIndicatorPoints(1));
      fprintf(fid,'%s,%d\n','ledY',registration_data.ledIndicatorPoints(2));
    end
  end
  
  fclose(fid);
  didsave = true;
catch ME
  warning('Could not save registration data to txt file: %s',getReport(ME));
end
if didsave,
  fprintf('Saved registration data to txt file %s\n',registrationtxtfile);
else
  fprintf('Could not save registration data to txt file:\n%s\n',getReport(ME));
end
end  % FlyDiscoRegisterTrx()



