function FlyDiscoDectectIndicatorLedOnOff(expdir,varargin)

% Deal with optional arguments
[analysis_protocol,settingsdir,datalocparamsfilestr] = ...
  myparse(varargin,...
  'analysis_protocol','current_bubble',...
  'settingsdir', default_settings_folder_path(), ...
  'datalocparamsfilestr','dataloc_params.txt');

% Write a header for this stage to the 'log'
timestamp = datestr(now,'yyyymmddTHHMMSS');
fprintf('\n\n***\nRunning %s with analysis protocol %s at %s\n', mfilename(), analysis_protocol, timestamp) ;

% Read in the data locations
analysis_protocol_folder_path = fullfile(settingsdir, analysis_protocol) ;
datalocparamsfile = fullfile(analysis_protocol_folder_path,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

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
registration_params = ReadParams(registrationparamsfile);
if ~isfield(registration_params,'OptogeneticExp') ,
  error('Registration params optogenetic flag %s does not exist',registrationparamsfile);
end
doLEDdetection = logical(registration_params.OptogeneticExp) ;
fprintf('\nOptogenetic flag set to: %g, running indicator detection\n',registration_params.OptogeneticExp);

% Load the registration data
registrationmatfile = fullfile(expdir,dataloc_params.registrationmatfilestr);
registration_data = load(registrationmatfile) ;

% Load metadata
metadata = collect_metadata(expdir, dataloc_params.metadatafilestr) ;

% Load trajectories
ctraxfile = fullfile(expdir,dataloc_params.ctraxfilestr);
moviefile = fullfile(expdir,dataloc_params.moviefilestr);
[trx,~,succeeded,timestamps] = load_tracks(ctraxfile,moviefile,'annname','');
if ~succeeded,
  error('Could not load trajectories from file %s',ctraxfile);
end

% Determine if timestamps are reliable, and compute a fallback dt if not
[are_timestamps_reliable, fallback_dt] = getFallbackDtIfNeeded(registration_params, timestamps, trx) ;

% Determine where the LED indicator is in the video frame
ledIndicatorPoints = ...
  detectLedLocations(registration_data, registration_params, ...
                     metadata, expdir, dataloc_params, timestamps, analysis_protocol_folder_path, ...
                     are_timestamps_reliable, fallback_dt) ;

% Make LED window, if an optogenetic experiment
if doLEDdetection ,
  indicatordata = extractIndicatorLED(expdir, dataloc_params, indicator_params, ledIndicatorPoints) ;
else
  indicatordata = struct() ;
end

% save indicator data to mat file
indicatordatamatfilename = fullfile(expdir,dataloc_params.indicatordatafilestr);
if exist(indicatordatamatfilename,'file'),
  delete(indicatordatamatfilename);
end
save(indicatordatamatfilename, '-struct', 'indicatordata') ;
fprintf('\n Saved indicator data to file %s\n', indicatordatamatfilename) ;

% Final message
fprintf('\n Finished running %s at %s.\n',mfilename(), datestr(now,'yyyymmddTHHMMSS'));
