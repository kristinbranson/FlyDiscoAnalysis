function FlyDiscoDetectIndicatorLedOnOff(expdir, varargin)

% Deal with optional arguments
[analysis_protocol, settingsdir, datalocparamsfilestr, forcecompute, debug, do_run] = ...
  myparse(varargin,...
  'analysis_protocol','current_bubble',...
  'settingsdir', default_settings_folder_path(), ...
  'datalocparamsfilestr','dataloc_params.txt', ...
  'forcecompute', false, ...
  'debug',false, ...
  'do_run', []) ;

% % Write a header for this stage to the 'log'
% timestamp = datestr(now,'yyyymmddTHHMMSS');
% fprintf('\n\n***\nRunning %s with analysis protocol %s at %s\n', mfilename(), analysis_protocol, timestamp) ;

% Read in the data locations
analysis_protocol_folder_path = fullfile(settingsdir, analysis_protocol) ;
datalocparamsfile = fullfile(analysis_protocol_folder_path, datalocparamsfilestr);
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
indicator_params = loadIndicatorParams(indicatorparamsfile) ;

% If the OptogeneticExp flag is false, issue a warning then exit, b/c they really
% shouldn't be running this stage in that case.  
if ~logical(indicator_params.OptogeneticExp) ,
  warning(['Optogenetic flag is set to false, but the ledonoffdetection stage is enabled.  ' ...
           'You should create an analysis-protocol-parameters.txt (if absent) in your analysis-protocol folder, and turn off this stage.']) ;
  % Save a "blank" indicator data to mat file, just to satisfy the stage's
  % required-files check.  We only do this for backwards-compatibility, and it
  % would be nice to remove at some point.
  indicatordata = struct() ;
  if exist(indicatordatamatfilename,'file'),
    delete(indicatordatamatfilename);
  end
  save(indicatordatamatfilename, '-struct', 'indicatordata') ;
  fprintf('Saved "empty" indicator data to file %s, just so the required-files check will pass\n', indicatordatamatfilename) ;  
  % Return early
  return
end

% Read in registration params
registrationparamsfile = fullfile(analysis_protocol_folder_path,dataloc_params.registrationparamsfilestr);
if ~exist(registrationparamsfile,'file')
  error('Registration params file %s does not exist',registrationparamsfile);
end
raw_registration_params = ReadParams(registrationparamsfile);
registration_params = modernizeRegistrationParams(raw_registration_params) ;

% Load the registration data
registrationmatfile = fullfile(expdir,dataloc_params.registrationmatfilestr);
registration_data = load(registrationmatfile) ;

% Load metadata
metadata = collect_metadata(expdir, dataloc_params.metadatafilestr) ;
rigId = metadata.rig ;  % Should be a scalar char array containing a single capital letter

% Load trajectories
ctraxfile = determine_downstream_trx_file_path(expdir, dataloc_params, do_run) ;
moviefile = fullfile(expdir,dataloc_params.moviefilestr);
[trx,~,succeeded,timestamps] = load_tracks(ctraxfile,moviefile,'annname','');
if ~succeeded,
  error('Could not load trajectories from file %s',ctraxfile);
end

% Get the frame interval for the movie
[are_timestamps_reliable, fallback_dt] = getFallbackDtIfNeeded(registration_params, timestamps, trx) ;
dt = fif(are_timestamps_reliable, median(diff(timestamps)), fallback_dt) ;  % frame interval, in seconds

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

% Make sure the protocol is shorter than the video
warn_if_protocol_is_longer_than_video(expdir, settingsdir, analysis_protocol, do_run) ;

% Determine where the LED indicator is in the video frame
[ledIndicatorPoints, templateShapeXY] = ...
  detectLedLocations(registration_data, ...
                     indicator_params, ...
                     expdir, ...
                     dataloc_params, ...
                     analysis_protocol_folder_path, ...
                     dt, ...
                     rigId) ;

% Save the LED detection image
if exist(imsavename,'file'),
  delete(imsavename);
end
saveLedDetectionImage(imsavename, bg_mean, ledIndicatorPoints, templateShapeXY, registration_params.doesYAxisPointUp) ;
fprintf('Saved LED detection image to file %s\n', imsavename) ;

% Read the LED masks, if present
has_explicit_masks = isfield(indicator_params, 'ledmask1') ;
if has_explicit_masks ,
  % Read in the masks
  for mask_index = 1:3 ,
    field_name = sprintf('ledmask%d', mask_index) ;
    if ~isfield(indicator_params, field_name) ,
      break
    end
    file_specification = indicator_params.(field_name) ;
    mask_file_path = determine_template_or_mask_file_path(file_specification, rigId, analysis_protocol_folder_path) ;
    raw_mask_as_double = double(imread(mask_file_path)) ;
    mask_as_double = mean(raw_mask_as_double, 3) ;
    threshold = 0.5*max(mask_as_double,[],'all') + 0.5*min(mask_as_double,[],'all') ;
    mask = (mask_as_double > threshold) ;
    mask_from_mask_index(:,:,mask_index) = mask ;  %#ok<AGROW> 
  end
else
  mask_from_mask_index = [] ;
end

% Extract indicatordata from the video
indicatordata = extractIndicatorLED(expdir, dataloc_params, indicator_params, ledIndicatorPoints, has_explicit_masks, mask_from_mask_index, debug) ;

% save indicator data to mat file
if exist(indicatordatamatfilename,'file'),
  delete(indicatordatamatfilename);
end
save(indicatordatamatfilename, '-struct', 'indicatordata') ;
fprintf('Saved indicator data to file %s\n', indicatordatamatfilename) ;

% Check if detected stim number matches number expected from protocol 
error_if_protocol_stim_num_notequal_detected(expdir, settingsdir, analysis_protocol, 'indicatordata', indicatordata); 

% % Final message
% fprintf('\n Finished running %s at %s.\n',mfilename(), datestr(now,'yyyymmddTHHMMSS'));

end
