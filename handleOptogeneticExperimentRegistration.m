function registration_data = ...
  handleOptogeneticExperimentRegistration(...
      registration_data, registration_params, metadata, expdir, dataloc_params, timestamps, analysis_protocol_folder_path)

% If not an optogenetic experiment, return immediately
if ~(isfield(registration_params,'OptogeneticExp') && registration_params.OptogeneticExp) ,
  return
end

% We want to modify the registration_params for handing to
% detectRegistrationMarks() for doing LED detection, so make a
% copy. (Not really needed, but it's nice to have the original values around
% for debugging.)
led_registration_params = registration_params ;

% At this point, we know led_registration_params.OptogeneticExp exists and that
% it's true


%
% For FlyBowlRGB and FlyBubbleRGB, convert RGB ledprotocol format to ChR led protocol format
%
didLoadProtocolOfSomeKind = false ;
if isfield(dataloc_params,'ledprotocolfilestr')
  led_protocol_file_path = fullfile(expdir,dataloc_params.ledprotocolfilestr) ;
  if exist(led_protocol_file_path,'file')
    protocolOfSomeKind = loadSingleVariableAnonymously(led_protocol_file_path, 'protocol') ;
    didLoadProtocolOfSomeKind = true ;
    if isExperimentRGB(metadata)
      if isfield(protocolOfSomeKind,'Rintensity')
        RGBprotocol = protocolOfSomeKind;
        % test if RGBprotocol has only one active color
        isLedActive = [any(RGBprotocol.Rintensity) any(RGBprotocol.Gintensity) any(RGBprotocol.Bintensity)] ;
        % check that there is 1 and only 1 color LED used in protocol
        activeLedCount = sum(isLedActive) ;
        if activeLedCount == 0
          error('ChR = 1 for LED protcol with no active LEDs')
        elseif activeLedCount > 1
          error('More than one active LED color in protocol. Not currently supported')
        end
        % call function that transforms new protocol to old protocol
        [protocol,ledcolor] = ConvertRGBprotocol2protocolformat(RGBprotocol, isLedActive);  %#ok<ASGLU>
      else
        protocol = protocolOfSomeKind ;
      end
    else
      protocol = protocolOfSomeKind ;
    end
  end
end


%
% Create max-value image for LED experiments  (needs fps)
%
if didLoadProtocolOfSomeKind ,
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
    if led_registration_params.usemediandt
      fps = 1/meddt;
    else
      fps = 1/median(diff(timestamps));
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
  if isfield(led_registration_params,'LEDMarkerType') && led_registration_params.OptogeneticExp,
    if ischar(led_registration_params.LEDMarkerType),
      if ~ismember(led_registration_params.LEDMarkerType,{'gradient'}),
        led_registration_params.LEDMarkerType = fullfile(analysis_protocol_folder_path,led_registration_params.LEDMarkerType);
      end
    else
      % rignames ABDC -> rigids
      % will need to change if it turns out to be cartridge
      % specific
      rigids = led_registration_params.LEDMarkerType(1:2:end-1);
      ledmarkertypes = led_registration_params.LEDMarkerType(2:2:end);
      rigid = metadata.rig;

      i = strcmp(rigid,rigids);
      if isempty(i),
        error('LEDMarkerType not set for plate %d',rigid);
      end
      if ~ismember(ledmarkertypes{i},{'gradient'}),
        led_registration_params.LEDMarkerType = fullfile(analysis_protocol_folder_path,ledmarkertypes{i});
      end
    end
  end

  % need to load find LEDMarkType for this plate
  if isfield(led_registration_params,'LEDMarkerType') && ischar(led_registration_params.LEDMarkerType) ,
    LEDimg = imread(led_registration_params.LEDMarkerType);
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


%
% Detect LED indicator using detectRegistrationMarks()
% 

% Tweak maxDistCornerFrac_LEDLabel, and stuff it into
% registration_params.maxDistCornerFrac_BowlLabel so that
% detectRegistrationMarks() can use it.
if isfield(led_registration_params,'maxDistCornerFrac_LEDLabel') ,
  led_registration_params.maxDistCornerFrac_BowlLabel = ...
    determineMaxDistCornerFracLabel(led_registration_params.maxDistCornerFrac_LEDLabel, metadata, 'maxDistCornerFrac_LEDLabel') ;
end

% Sub in the LED marker type fpr the bowlMarkerType
led_registration_params.bowlMarkerType = led_registration_params.LEDMarkerType;

% Collect the registration params into a cell array for passing to
% detectRegistrationMarks()
led_registration_params_cell = ...
  marshallRegistrationParamsForDetectRegistrationMarks(led_registration_params, expdir, dataloc_params.ledregistrationimagefilstr) ;

% detect
ledindicator_data = ...
  detectRegistrationMarks(...
    led_registration_params_cell{:}, ...
    'bkgdImage', im2double(ledMaxImage), ...
    'ledindicator', true, ...
    'regXY', registration_data.bowlMarkerPoints, ...
    'useNormXCorr',true);
registration_data.ledIndicatorPoints = ledindicator_data.bowlMarkerPoints;

% Decare victory
fprintf('Detected led indicator.\n');

end  % function

