function [success,msgs] = FlyDiscoRegisterTrx(expdir,varargin)

success = true;
msgs = {};

%version = '0.2';

fns_notperframe = {'id','moviename','annname','firstframe','arena','off',...
  'nframes','endframe','matname','fps','pxpermm'};

[analysis_protocol,settingsdir,registrationparamsfilestr,datalocparamsfilestr,dotemporalreg] = ...
  myparse(varargin,...
  'analysis_protocol','current_bubble',...
  'settingsdir','/groups/branson/home/robiea/Code_versioned/FlyBubbleAnalysis/settings',...
  'registrationparamsfilestr','registration_params.txt',...
  'datalocparamsfilestr','dataloc_params.txt',...
  'dotemporalreg',false);

%% read in the data locations
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

% logger = PipelineLogger(expdir,mfilename(),dataloc_params,...
%   'registertrx_logfilestr',settingsdir,analysis_protocol,...
%   'versionstr',version);

%% read in registration params

% name of parameters file
registrationparamsfile = fullfile(settingsdir,analysis_protocol,registrationparamsfilestr);
if ~exist(registrationparamsfile,'file'),
  error('Registration params file %s does not exist',registrationparamsfile);
end
% read
registration_params = ReadParams(registrationparamsfile);
if isfield(registration_params,'doTemporalRegistration'),
  dotemporalreg = registration_params.doTemporalRegistration;
end

%% detect registration marks

% name of annotation file
%annfile = fullfile(expdir,dataloc_params.annfilestr);
annfile = '';
bg_mean = [];
if isfield(dataloc_params,'flytrackerbgstr'),
  try
    flytrackerbgfile = fullfile(expdir,dataloc_params.flytrackerbgstr);
    load(flytrackerbgfile,'bg');
    bg_mean = 255*bg.bg_mean;
    clear bg;
  catch ME,
    warning('Could not load background image from %s',flytrackerbgfile);
    disp(getReport(ME));
  end
end
    

% name of movie file
moviefile = fullfile(expdir,dataloc_params.moviefilestr);

% Collect metadata
metadata = collect_metadata(expdir, dataloc_params.metadatafilestr) ;

% template filename should be relative to settings directory
if isfield(registration_params,'bowlMarkerType'),
  if ischar(registration_params.bowlMarkerType),
    if ~ismember(registration_params.bowlMarkerType,{'gradient'}),
      registration_params.bowlMarkerType = fullfile(settingsdir,analysis_protocol,registration_params.bowlMarkerType);
    end
  else
    % plate -> bowlmarkertype
    plateids = registration_params.bowlMarkerType(1:2:end-1);
    bowlmarkertypes = registration_params.bowlMarkerType(2:2:end);
%     [metadata,success1] = parseExpDir(expdir);
%     if ~success1 || ~isfield(metadata,'plate'),
%       metadata = ReadMetadataFile(fullfile(expdir,dataloc_params.metadatafilestr));
%     end
    if isnumeric(metadata.plate),
        plateid = num2str(metadata.plate);
    else
        plateid = metadata.plate;
    end
    if iscell(plateids)
        i = find(strcmp(num2str(plateid), plateids));
    else
        i = find(str2double(plateid) == plateids,1);
    end
    if isempty(i),
      error('bowlMarkerType not set for plate %s',plateid);
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
%   [metadata,success1] = parseExpDir(expdir);
%   if ~success1 || ~isfield(metadata,'plate'),
%     metadata = ReadMetadataFile(fullfile(expdir,dataloc_params.metadatafilestr));
%   end
  if isnumeric(metadata.plate),
    plateid = num2str(metadata.plate);
  else
    plateid = metadata.plate;
  end
  if iscell(plateids)
      i = find(strcmp(num2str(plateid), plateids));
  else
      i = find(str2double(plateid) == plateids,1);
  end
  if isempty(i),
    error('maxDistCornerFrac_BowlLabel not set for plate %d',plateid);
  end
  if iscell(cornerfracs)
      registration_params.maxDistCornerFrac_BowlLabel = str2double(cornerfracs{i});
  else
      registration_params.maxDistCornerFrac_BowlLabel = cornerfracs(i);
  end
end

fnsignore = intersect(fieldnames(registration_params),...
  {'minFliesLoadedTime','maxFliesLoadedTime','extraBufferFliesLoadedTime','usemediandt','doTemporalRegistration','OptogeneticExp','LEDMarkerType','maxDistCornerFrac_LEDLabel'});
  
registration_params_cell = struct2paramscell(rmfield(registration_params,fnsignore));
% file to save image to
if isfield(dataloc_params,'registrationimagefilestr'),
  registration_params_cell(end+1:end+2) = {'imsavename',fullfile(expdir,dataloc_params.registrationimagefilestr)};
end
% detect
try
  registration_data = detectRegistrationMarks(registration_params_cell{:},'bkgdImage',bg_mean,'movieName',moviefile);
catch ME,
  fprintf('Error detecting registration marks:\n');
  fprintf(getReport(ME));
  success = false;
  msgs = {['Error detecting registration marks: ',getReport(ME)]};
  return;
end

fprintf('Detected registration marks.\n');

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

% remove nans from flytracker outputs
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

registration_data.flytracker_nnanframes = ninterpframes;
registration_data.flytracker_nids0 = nids0;

% frame rate
if registration_params.usemediandt,
  tmp = diff(timestamps);
  tmp(isnan(tmp)) = [];
  if isempty(tmp),
    meddt = 1/trx(1).fps;
  else
    meddt = median(tmp);
  end
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
  if registration_params.usemediandt,
    trx(fly).dt = repmat(meddt,[1,trx(fly).nframes-1]);
  else
    if isfield(trx,'timestamps'),
      trx(fly).dt = diff(trx(fly).timestamps);
    else
      trx(fly).dt = repmat(1/trx(fly).fps,[1,trx(fly).nframes-1]);
    end
  end

  trx(fly).fps = 1/mean(trx(fly).dt);
  trx(fly).pxpermm = 1 / registration_data.scale;
  
end

fprintf('Applied spatial registration.\n');
%% for FlyBowlRGB and FlyBubbleRGB convert RGB ledprotocol format to ChR led protocol format
if isfield(registration_params,'OptogeneticExp')
    if registration_params.OptogeneticExp
        if isfield(dataloc_params,'ledprotocolfilestr')
            if exist(fullfile(expdir,dataloc_params.ledprotocolfilestr),'file')
                load(fullfile(expdir,dataloc_params.ledprotocolfilestr),'protocol')
                if strcmp(metadata.assay,'FlyBubbleRGB') || strcmp(metadata.assay,'FlyBowlRGB')                    
                    if isfield(protocol,'Rintensity')
                        RGBprotocol = protocol;
                        clear protocol;
                        % test if RGBprotocol has only one active color
                        countactiveLEDs = [double(any(RGBprotocol.Rintensity));double(any(RGBprotocol.Gintensity));double(any(RGBprotocol.Bintensity))];
                        % check that there is 1 and only 1 color LED used in protocol
                        if sum(countactiveLEDs) == 0
                            error('ChR = 1 for LED protcol with no active LEDs')
                        elseif sum(countactiveLEDs) > 1
                            error('More than one active LED color in protocol. Not currently supported')
                        end
                        % call function that transforms new protocol to old protocol
                        [protocol,ledcolor] = ConvertRGBprotocol2protocolformat(RGBprotocol,countactiveLEDs);
                    end
                end
            end
        end
    end
end

%% create maxvalueimage for LED experiments  (needs fps)
if isfield(registration_params,'OptogeneticExp') ,
    if registration_params.OptogeneticExp && exist('protocol','var')
        moviefilename = fullfile(expdir,dataloc_params.moviefilestr);
        [readfcn,~,~,headerinfo] = get_readframe_fcn(moviefilename);
        %determine minimum frames to process
        % protocol now loaded in earlier cell - converted if need for RGB
        firstactiveStepNum = find([protocol.intensity] ~= 0,1);
        if isempty(firstactiveStepNum)
            error('ChR = 1 for LED protcol with no active LEDs')
        end
        
        % experiments
        if ~isempty(protocol)
            fps = 1/meddt;
            secstoLEDpulse = protocol.delayTime(firstactiveStepNum) + protocol.pulsePeriodSP(firstactiveStepNum)*protocol.pulseNum(firstactiveStepNum)/1000;
            % sets end range in which to find LED on
            frametoLEDpulse = secstoLEDpulse*fps;
            if protocol.pulsePeriodSP(firstactiveStepNum)*protocol.pulseNum(firstactiveStepNum) <= 1/fps*1000*2 %pulseWidth(ms) < 2 times sampling
                jump = 1;
            else
                % in frames
                jump = round(protocol.pulsePeriodSP(firstactiveStepNum)*protocol.pulseNum(firstactiveStepNum)/1000*fps/2);
            end
            
        else
            % reasonable guess
            frametoLEDpulse = headerinfo.nframes/3;
            jump = 10;
        end
        frametoLEDpulse = max([round(frametoLEDpulse),min(200,headerinfo.nframes)]);
        jump = min(jump,round(frametoLEDpulse/100));
        im = readfcn(1);
        for j = 1:jump:frametoLEDpulse
            tmp = readfcn(j);
            idx = tmp > im;
            im(idx) = tmp(idx);
            %   if (mod(j,1000) == 0);
            %     disp(round((j/frametoLEDpulse)*100))
            %   end
        end
        % check if im has indicator on
        if isfield(registration_params,'LEDMarkerType') && registration_params.OptogeneticExp,
            if ischar(registration_params.LEDMarkerType),
                if ~ismember(registration_params.LEDMarkerType,{'gradient'}),
                    registration_params.LEDMarkerType = fullfile(settingsdir,analysis_protocol,registration_params.LEDMarkerType);
                end
            else
                % rignames ABDC -> rigids
                % will need to change if it turns out to be cartridge
                % specific
                rigids = registration_params.LEDMarkerType(1:2:end-1);
                ledmarkertypes = registration_params.LEDMarkerType(2:2:end);
%                 [metadata,success1] = parseExpDir(expdir);
%                 if ~success1 || ~isfield(metadata,'rig'),
%                     metadata = ReadMetadataFile(fullfile(expdir,dataloc_params.metadatafilestr));
%                 end
                rigid = metadata.rig;
                
                i = strcmp(rigid,rigids);
                if isempty(i),
                    error('LEDMarkerType not set for plate %d',rigid);
                end
                if ~ismember(ledmarkertypes{i},{'gradient'}),
                    registration_params.LEDMarkerType = fullfile(settingsdir,analysis_protocol,ledmarkertypes{i});
                end
            end
        end
        
        % need to load find LEDMarkType for this plate
        if isfield(registration_params,'LEDMarkerType') && ischar(registration_params.LEDMarkerType) && registration_params.OptogeneticExp,
            LEDimg = imread(registration_params.LEDMarkerType);
            binLEDstep = 25;
            hist_LED = histcounts(LEDimg,(0:binLEDstep:255));
            brightthres = sum(hist_LED(9:10));
            diffimage = im - readfcn(1);
            
            [xgrid, ygrid] = meshgrid(1:size(diffimage,2), 1:size(diffimage,1));
            mask = ((xgrid-registration_data.circleCenterX).^2 + (ygrid-registration_data.circleCenterY).^2) >= registration_data.circleRadius.^2;
            outSideArenaPxs = diffimage(mask);
            hist_outSideArenaPxs = histcounts(outSideArenaPxs,(0:binLEDstep:255));
            brightpxs = sum(hist_outSideArenaPxs(9:10));
            if brightpxs <= brightthres %|| brightpxs >= 15
                % no LED indicator in im, redo for whole movie
                im = readfcn(1);
                for j = 1:jump:headerinfo.nframes
                    tmp = readfcn(j);
                    idx = tmp > im;
                    im(idx) = tmp(idx);
                end
            end
        end
    end
end
%% detect LED indicator (modified from detect registration marks)
if isfield(registration_params,'OptogeneticExp')
    if registration_params.OptogeneticExp
        % name of annotation file
        annfile = fullfile(expdir,dataloc_params.annfilestr);
        
        % name of movie file
        moviefile = fullfile(expdir,dataloc_params.moviefilestr);
        
        % LED MarkerType set to metadata specific rig in previous cell 
        
%         if isfield(registration_params,'LEDMarkerType'),
%             if ischar(registration_params.LEDMarkerType),
%                 if ~ismember(registration_params.LEDMarkerType,{'gradient'}),
%                     registration_params.LEDMarkerType = fullfile(settingsdir,analysis_protocol,registration_params.LEDMarkerType);
%                 end
%             else
%                 % rignames ABDC -> rigids
%                 % will need to change if it turns out to be cartridge
%                 % specific
%                 rigids = registration_params.LEDMarkerType(1:2:end-1);
%                 ledmarkertypes = registration_params.LEDMarkerType(2:2:end);
%                 [metadata,success1] = parseExpDir(expdir);
%                 if ~success1 || ~isfield(metadata,'rig'),
%                     metadata = ReadMetadataFile(fullfile(expdir,dataloc_params.metadatafilestr));
%                 end
%                 rigid = metadata.rig;
%                 
%                 i = strcmp(rigid,rigids);
%                 if isempty(i),
%                     error('LEDMarkerType not set for plate %d',rigid);
%                 end
%                 if ~ismember(ledmarkertypes{i},{'gradient'}),
%                     registration_params.LEDMarkerType = fullfile(settingsdir,analysis_protocol,ledmarkertypes{i});
%                 end
%             end
%         end
        
        % maxDistCornerFrac_BowlLabel might depend on bowl
        if isfield(registration_params,'maxDistCornerFrac_LEDLabel') && ...
            numel(registration_params.maxDistCornerFrac_LEDLabel) > 1,
            plateids = registration_params.maxDistCornerFrac_LEDLabel(1:2:end-1);
            cornerfracs = registration_params.maxDistCornerFrac_LEDLabel(2:2:end);
%             [metadata,success1] = parseExpDir(expdir);
%             if ~success1 || ~isfield(metadata,'plate'),
%                 metadata = ReadMetadataFile(fullfile(expdir,dataloc_params.metadatafilestr));
%             end
            if isnumeric(metadata.plate),
                plateid = num2str(metadata.plate);
            else
                plateid = metadata.plate;
            end
            if iscell(plateids)
                i = find(strcmp(num2str(plateid), plateids));
            else
                i = find(str2double(plateid) == plateids,1);
            end
            if isempty(i),
                error('maxDistCornerFrac_LEDLabel not set for plate %d',plateid);
            end
            if iscell(cornerfracs)
                registration_params.maxDistCornerFrac_BowlLabel = str2double(cornerfracs{i});
            else
                registration_params.maxDistCornerFrac_BowlLabel = cornerfracs(i);
            end
        else
            
            registration_params.maxDistCornerFrac_BowlLabel = registration_params.maxDistCornerFrac_LEDLabel;
            
        end
        
        
        
        registration_params.bowlMarkerType = registration_params.LEDMarkerType;
       
        
        fnsignore = intersect(fieldnames(registration_params),...
            {'minFliesLoadedTime','maxFliesLoadedTime','extraBufferFliesLoadedTime','usemediandt','doTemporalRegistration','OptogeneticExp','LEDMarkerType','maxDistCornerFrac_LEDLabel'});
        
        registration_params_cell = struct2paramscell(rmfield(registration_params,fnsignore));
        
       
        % file to save image to
        if isfield(dataloc_params,'ledregistrationimagefilstr'),
            registration_params_cell(end+1:end+2) = {'imsavename',fullfile(expdir,dataloc_params.ledregistrationimagefilstr)};
        end
        
        
        % detect
        try
            ledindicator_data = detectRegistrationMarks(registration_params_cell{:},'bkgdImage',bg_mean,'movieName',moviefile,'bkgdImage',im2double(im),'ledindicator',true,'regXY',registration_data.bowlMarkerPoints);
            
        catch ME,
            fprintf('Error detecting led indicator:\n');
            fprintf(getReport(ME));
            success = false;
            msgs = {['Error detecting led indicator: ',getReport(ME)]};
            return;
        end
        registration_data.ledIndicatorPoints = ledindicator_data.bowlMarkerPoints;
        
        fprintf('Detected led indicator.\n');
               
        
    end
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
    newid2oldid(trxdelete) = [];
    
  end
  
  fprintf('Applied temporal registration.\n');
  
else
  
  fprintf('NOT applying temporal registration.\n');
  
end

%% save registered trx to file

% name of output trx mat file
trxfile = fullfile(expdir,dataloc_params.trxfilestr);

didsave = false;
try
  if exist(trxfile,'file'),
    delete(trxfile);
  end
  %registrationinfo = logger.runInfo;
  save(trxfile,'trx','timestamps');
  didsave = true;
catch ME
  warning('Could not save registered trx: %s',getReport(ME));
  success = false;
  msgs{end+1} = ['Could not save registered trx: %s',getReport(ME)];
end

if didsave,
  fprintf('Saved registered trx to file %s\n',trxfile);
else
  fprintf('Could not save registered trx:\n%s\n',getReport(ME));
end

%% save params to mat file

registration_data.newid2oldid = newid2oldid;
registration_data.flytracker_nidsnew = numel(trx);

registrationmatfile = fullfile(expdir,dataloc_params.registrationmatfilestr);
tmp = rmfield(registration_data,'registerfn'); 
%tmp.registrationinfo = registrationinfo;

didsave = false;
try
  if exist(registrationmatfile,'file'),
    delete(registrationmatfile);
  end
  save(registrationmatfile,'-struct','tmp');
  didsave = true;
catch ME
  warning('Could not save registered data to mat file: %s',getReport(ME));
  success = false;
  msgs{end+1} = ['Could not save registered data to mat file: %s',getReport(ME)];
end

if didsave,
  fprintf('Saved registration data to file %s\n',registrationmatfile);
else
  fprintf('Could not save registration data to mat file:\n%s\n',getReport(ME));
end

%% save params to text file
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
%   fnssave = {'version','timestamp','analysis_protocol','linked_analysis_protocol'};
%   fnssave = intersect(fnssave,fieldnames(registrationinfo));
%   for i = 1:numel(fnssave),
%     fn = fnssave{i};
%     fprintf(fid,'%s,%s\n',fn,registrationinfo.(fn));
%   end
  
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
  success = false;
  msgs{end+1} = ['Could not save registration data to txt file: %s',getReport(ME)];
end

if didsave,
  fprintf('Saved registration data to txt file %s\n',registrationtxtfile);
else
  fprintf('Could not save registration data to txt file:\n%s\n',getReport(ME));
end

% %%
% logger.close();

%% close figures

if isdeployed,
  delete(findall(0,'type','figure'));
end
