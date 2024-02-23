function result = ...
  detectLedLocations(...
      registration_data, registration_params, expdir, dataloc_params, analysis_protocol_folder_path, ...
      bkgdImage, ...
      dt, ...
      rigId)

% If not an optogenetic experiment, error
if ~(isfield(registration_params,'OptogeneticExp') && registration_params.OptogeneticExp) ,
  error('LED detection requested on a non-OptogeneticExp experiment') ;
end

% At this point, we know registration_params.OptogeneticExp exists and that
% it's true


%
% For FlyBowlRGB and FlyBubbleRGB, convert RGB ledprotocol format to ChR led protocol format
%
if isfield(dataloc_params,'ledprotocolfilestr')
  led_protocol_file_path = fullfile(expdir,dataloc_params.ledprotocolfilestr) ;
else
  error('No ledprotocolfilestr defined in dataloc_params') ;
end
if exist(led_protocol_file_path,'file')
  protocolOfSomeKind = loadSingleVariableAnonymously(led_protocol_file_path, 'protocol') ;
  protocol = downmixProtocolIfNeeded(protocolOfSomeKind) ;
else
  error('LED protocol file %s does not exist', dataloc_params.ledprotocolfilestr) ;
end


%
% Create max-value image for LED experiments, and get the LED template image
%
[ledMaxImage, template] = ...
  determineLedMaxImageAndLedTemplate(expdir, analysis_protocol_folder_path, dataloc_params, registration_params, registration_data, protocol, ...
                                     dt, rigId) ;



%
% Detect LED indicator using findTemplateMatchWithPossibleRotation()
% 

% Find the location of the LED(s)
ledMaxImageDouble = im2double(ledMaxImage) ;
minTemplateFeatureStrength = registration_params.minTemplateFeatureStrength ;
nRotations = 20 ;
useNormXCorr = true ;
[nr,nc,~] = size(ledMaxImageDouble) ;
excluded_radius = registration_params.maxDistCornerFrac_LEDLabel * min(nr,nc) ;
excluded_xyrs = [registration_data.bowlMarkerPoints(1) registration_data.bowlMarkerPoints(2) excluded_radius] ;
corner_radius = registration_params.maxDistCornerFrac_LEDLabel * min(nr,nc) ;
result = findTemplateMatchWithPossibleRotation(ledMaxImageDouble, ...
                                               template, ...
                                               minTemplateFeatureStrength, ...
                                               nRotations, ...
                                               useNormXCorr, ...
                                               corner_radius, ...
                                               excluded_xyrs) ;

% Decare victory
fprintf('Detected led indicator.\n');

% Save the LED detection image
imsavename = fullfile(expdir, dataloc_params.ledregistrationimagefilstr) ;
saveLedDetectionImage(imsavename, bkgdImage, result) ;

end  % function

