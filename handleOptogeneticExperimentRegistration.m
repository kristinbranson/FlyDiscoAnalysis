function registration_params = handleOptogeneticExperimentRegistration(registration_params, expdir, dataloc_params)

% If not an optogenetic experiment, return immediately
if ~(isfield(registration_params,'OptogeneticExp') && registration_params.OptogeneticExp) ,
  return
end

% At this point, we know registration_params.OptogeneticExp exists and that
% it's true

%%% For FlyBowlRGB and FlyBubbleRGB convert RGB ledprotocol format to ChR led protocol format
if isfield(dataloc_params,'ledprotocolfilestr')
  led_protocol_file_path = fullfile(expdir,dataloc_params.ledprotocolfilestr) ;
  if exist(led_protocol_file_path,'file')
    protocolOfSomeKind = loadSingleVariableAnonymously(led_protocol_file_path, 'protocol') ;
    if isExperimentRGB(metadata)
      if isfield(protocolOfSomeKind,'Rintensity')
        RGBprotocol = protocolOfSomeKind;
        % test if RGBprotocol has only one active color
        countactiveLEDs = [double(any(RGBprotocol.Rintensity));double(any(RGBprotocol.Gintensity));double(any(RGBprotocol.Bintensity))];
        % check that there is 1 and only 1 color LED used in protocol
        if sum(countactiveLEDs) == 0
          error('ChR = 1 for LED protcol with no active LEDs')
        elseif sum(countactiveLEDs) > 1
          error('More than one active LED color in protocol. Not currently supported')
        end
        % call function that transforms new protocol to old protocol
        [protocol,ledcolor] = ConvertRGBprotocol2protocolformat(RGBprotocol,countactiveLEDs);  %#ok<ASGLU>
      else
        protocol = protocolOfSomeKind ;
      end
    else
      protocol = protocolOfSomeKind ;
    end
  end
end

%%% Create max-value image for LED experiments  (needs fps)
if exist('protocol','var')
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
    if registration_params.usemediandt
      fps = 1/meddt;
    else
      if isvarname('timestamps')
        fps = 1/median(diff(timestamps));
      else
        [~,~,succeeded,timestamps] = load_tracks(ctraxfile,moviefile);
        if ~succeeded
          error('Could not load trajectories from file %s',ctraxfile);
        end
        fps = 1/median(diff(timestamps));
      end
    end

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
  ledMaxImage = sampleFramesForMaximumImage(readfcn, 1, jump, frametoLEDpulse) ;
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
    diffimage = ledMaxImage - readfcn(1);

    [xgrid, ygrid] = meshgrid(1:size(diffimage,2), 1:size(diffimage,1));
    mask = ((xgrid-registration_data.circleCenterX).^2 + (ygrid-registration_data.circleCenterY).^2) >= registration_data.circleRadius.^2;
    outSideArenaPxs = diffimage(mask);
    hist_outSideArenaPxs = histcounts(outSideArenaPxs,(0:binLEDstep:255));
    brightpxs = sum(hist_outSideArenaPxs(9:10));
    if brightpxs <= brightthres
      % no LED indicator in im, redo for whole movie
      ledMaxImage = sampleFramesForMaximumImage(readfcn, 1, jump, headerinfo.nframes) ;
    end
  end
end

%%% detect LED indicator (modified from detect registration marks)
% name of movie file
moviefile = fullfile(expdir,dataloc_params.moviefilestr);

% maxDistCornerFrac_BowlLabel might depend on bowl
if isfield(registration_params,'maxDistCornerFrac_LEDLabel') && ...
    numel(registration_params.maxDistCornerFrac_LEDLabel) > 1,
  plateids = registration_params.maxDistCornerFrac_LEDLabel(1:2:end-1);
  cornerfracs = registration_params.maxDistCornerFrac_LEDLabel(2:2:end);
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

fnsignore = ...
  intersect(fieldnames(registration_params),...
  {'minFliesLoadedTime', ...
  'maxFliesLoadedTime', ...
  'extraBufferFliesLoadedTime', ...
  'usemediandt', ...
  'doTemporalRegistration', ...
  'OptogeneticExp', ...
  'LEDMarkerType', ...
  'maxDistCornerFrac_LEDLabel', ...
  'doTemporalTruncation', ...
  'maxFlyTrackerNanInterpFrames', ...
  'bkgdNSampleFrames'});
registration_params_cell = struct2paramscell(rmfield(registration_params,fnsignore));

% file to save image to
if isfield(dataloc_params,'ledregistrationimagefilstr'),
  registration_params_cell(end+1:end+2) = {'imsavename',fullfile(expdir,dataloc_params.ledregistrationimagefilstr)};
end

% detect
ledindicator_data = ...
  detectRegistrationMarks(registration_params_cell{:}, ...
  'bkgdImage', im2double(ledMaxImage), ...
  'ledindicator', true, ...
  'regXY', registration_data.bowlMarkerPoints, ...
  'useNormXCorr',true);
registration_data.ledIndicatorPoints = ledindicator_data.bowlMarkerPoints;

% Decare victory
fprintf('Detected led indicator.\n');

end  % function

