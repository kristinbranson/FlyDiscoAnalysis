function FlyDiscoDetectIndicatorLedOnOff(expdir, varargin)

% Deal with optional arguments
[analysis_protocol, settingsdir, datalocparamsfilestr, forcecompute, debug] = ...
  myparse(varargin,...
  'analysis_protocol','current_bubble',...
  'settingsdir', default_settings_folder_path(), ...
  'datalocparamsfilestr','dataloc_params.txt', ...
  'forcecompute', false, ...
  'debug',false) ;

% Write a header for this stage to the 'log'
timestamp = datestr(now,'yyyymmddTHHMMSS');
fprintf('\n\n***\nRunning %s with analysis protocol %s at %s\n', mfilename(), analysis_protocol, timestamp) ;

% Read in the data locations
analysis_protocol_folder_path = fullfile(settingsdir, analysis_protocol) ;
datalocparamsfile = fullfile(analysis_protocol_folder_path,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

% If forcecompute is true, delete the output files now
imsavename = fullfile(expdir, dataloc_params.ledregistrationimagefilstr) ;
indicatordatamatfilename = fullfile(expdir,dataloc_params.indicatordatafilestr);
if forcecompute ,
  % Delete the output files, so they have to be recomputed.
  % This may not particularly matter for LED indicator detection, actually...

  % Delete the LED detection image
  if exist(imsavename,'file'),
    delete(imsavename);
  end

  % Delete indicator data mat file
  if exist(indicatordatamatfilename,'file'),
    delete(indicatordatamatfilename);
  end
end

% Read in indicator params
indicatorparamsfile = fullfile(analysis_protocol_folder_path,dataloc_params.indicatorparamsfilestr);
if ~exist(indicatorparamsfile,'file')
  error('Indictor params file %s does not exist',indicatorparamsfile);
end
indicator_params = ReadParams(indicatorparamsfile);

% Read in registration params and check for optogenetic flag
registrationparamsfile = fullfile(analysis_protocol_folder_path,dataloc_params.registrationparamsfilestr);
if ~exist(registrationparamsfile,'file')
  error('Registration params file %s does not exist',registrationparamsfile);
end
raw_registration_params = ReadParams(registrationparamsfile);
registration_params = modernizeRegistrationParams(raw_registration_params) ;
if ~isfield(indicator_params,'OptogeneticExp') ,
  error('Indicator params optogenetic flag %s does not exist',registrationparamsfile);
end
doLEDdetection = logical(indicator_params.OptogeneticExp) ;
if doLEDdetection
  fprintf('\nOptogenetic flag set to: %g, running indicator detection\n',indicator_params.OptogeneticExp);
else
  fprintf('\nOptogenetic flag set to: %g, not running indicator detection\n',indicator_params.OptogeneticExp);
end

% Load the registration data
registrationmatfile = fullfile(expdir,dataloc_params.registrationmatfilestr);
registration_data = load(registrationmatfile) ;

% Load metadata
metadata = collect_metadata(expdir, dataloc_params.metadatafilestr) ;
rigId = metadata.rig ;  % Should be a scalar char array containing a single capital letter

% Load trajectories
ctraxfile = fullfile(expdir,dataloc_params.ctraxfilestr);
moviefile = fullfile(expdir,dataloc_params.moviefilestr);
[trx,~,succeeded,timestamps] = load_tracks(ctraxfile,moviefile,'annname','');
if ~succeeded,
  error('Could not load trajectories from file %s',ctraxfile);
end

% Get the frame interval for the movie
[are_timestamps_reliable, fallback_dt] = getFallbackDtIfNeeded(registration_params, timestamps, trx) ;
dt = fif(are_timestamps_reliable, median(diff(timestamps)), fallback_dt) ;  % frame interval, in seconds

% Make LED window, if an optogenetic experiment
if doLEDdetection ,
  % Load the background file output by FlyTracker
  if isfield(dataloc_params,'flytrackerbgstr')
    flytrackerbgfile = fullfile(expdir,dataloc_params.flytrackerbgstr);
    load(flytrackerbgfile,'bg');
    bg_mean = 255*bg.bg_mean;  % a double array, elements on [0,255], not necessarily integers
  else
    % Note: We really do want to error if this is missing.
    % We've been bitten by this not being what we thought it was.
    error('dataloc_params is missing field flytrackerbgstr, which is required');
  end
  
  % Determine where the LED indicator is in the video frame
  ledIndicatorPoints = ...
    detectLedLocations(registration_data, ...
                       indicator_params, ...
                       expdir, ...
                       dataloc_params, ...
                       analysis_protocol_folder_path, ...
                       dt, ...
                       rigId) ;
  % Extract indicatordata from the video
  indicatordata = extractIndicatorLED(expdir, dataloc_params, indicator_params, ledIndicatorPoints, debug) ;
else
  indicatordata = struct() ;
end

% Save the LED detection image
if exist(imsavename,'file'),
  delete(imsavename);
end
saveLedDetectionImage(imsavename, bg_mean, ledIndicatorPoints, registration_params.doesYAxisPointUp) ;
fprintf('Saved LED detection image to file %s\n', imsavename) ;

% save indicator data to mat file
if exist(indicatordatamatfilename,'file'),
  delete(indicatordatamatfilename);
end
save(indicatordatamatfilename, '-struct', 'indicatordata') ;
fprintf('Saved indicator data to file %s\n', indicatordatamatfilename) ;

% Final message
fprintf('\n Finished running %s at %s.\n',mfilename(), datestr(now,'yyyymmddTHHMMSS'));
