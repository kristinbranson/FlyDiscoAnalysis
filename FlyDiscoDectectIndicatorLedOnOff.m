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
datalocparamsfile = fullfile(settingsdir,analysis_protocol,datalocparamsfilestr);
dataloc_params = ReadParams(datalocparamsfile);

% Read in indicator params
indicatorparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.indicatorparamsfilestr);
if ~exist(indicatorparamsfile,'file')
  error('Indictor params file %s does not exist',indicatorparamsfile);
end
indicator_params = ReadParams(indicatorparamsfile);

% Read in registration params and check for optogenetic flag
registrationparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.registrationparamsfilestr);
if ~exist(registrationparamsfile,'file')
  error('Registration params file %s does not exist',registrationparamsfile);
end
registration_params = ReadParams(registrationparamsfile);
if ~isfield(registration_params,'OptogeneticExp') ,
  error('Registration params optogenetic flag %s does not exist',registrationparamsfile);
end
doLEDdetection = logical(registration_params.OptogeneticExp) ;
fprintf('\nOptogenetic flag set to: %g, running indicator detection\n',registration_params.OptogeneticExp);

% Make LED window, if an optogenetic experiment
if doLEDdetection ,
  indicatordata = extractIndicatorLED(expdir, dataloc_params, indicator_params) ;
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
