function [result, templateShapeXY] = ...
  detectLedLocations(...
      registration_data, ...
      indicator_params, ...
      expdir, ...
      dataloc_params, ...
      analysis_protocol_folder_path, ...
      dt, ...
      rigId)

% This should only be used with optogenetic experiments

%
% For FlyBowlRGB and FlyBubbleRGB, convert RGB ledprotocol format to ChR led protocol format
%
if indicator_params.doesUseProtocolDotMat
  if isfield(dataloc_params,'ledprotocolfilestr')
    led_protocol_file_path = fullfile(expdir, dataloc_params.ledprotocolfilestr) ;
  else
    error('No ledprotocolfilestr defined in dataloc_params') ;
  end
  if exist(led_protocol_file_path,'file')
    protocolOfSomeKind = loadSingleVariableAnonymously(led_protocol_file_path, 'protocol') ;
    protocolOrEmpty = downmixProtocolIfNeeded(protocolOfSomeKind, indicator_params) ;
  else
    warning('LED protocol file %s does not exist, so using best guess about frames in which to look for lighted LED', dataloc_params.ledprotocolfilestr) ;
    protocolOrEmpty = [] ;
  end
else
  protocolOrEmpty = [] ;
end  

%
% Create max-value image for LED experiments, and get the LED template image
%
[ledMaxImage, template] = ...
  determineLedMaxImageAndLedTemplate(expdir, ...
                                     analysis_protocol_folder_path, ...
                                     dataloc_params, ...
                                     indicator_params, ...
                                     registration_data, ...
                                     protocolOrEmpty, ...
                                     dt, ...
                                     rigId) ;
templateShapeXY = fliplr(size(template)) ;


%
% Detect LED indicator using findTemplateMatchWithPossibleRotation()
% 

% Find the location of the LED(s)
ledMaxImageDouble = im2double(ledMaxImage) ;
minTemplateFeatureStrength = -0.5 ;
if isfield(indicator_params, 'nRotations')
  nRotations = indicator_params.nRotations ;
else
  nRotations = 20 ;
end
if isfield(indicator_params, 'useNormXCorr')
  useNormXCorr = indicator_params.useNormXCorr ;
else
  useNormXCorr = true ;
end
[nr,nc,~] = size(ledMaxImageDouble) ;
excluded_radius = indicator_params.maxDistCornerFrac_LEDLabel * min(nr,nc) ;
if isempty(registration_data.bowlMarkerPoints)
  excluded_xyrs = zeros(0,3) ;
else
  excluded_xyrs = [registration_data.bowlMarkerPoints(1) registration_data.bowlMarkerPoints(2) excluded_radius] ;
end
corner_radius = indicator_params.maxDistCornerFrac_LEDLabel * min(nr,nc) ;
result = findTemplateMatchWithPossibleRotation(ledMaxImageDouble, ...
                                               template, ...
                                               minTemplateFeatureStrength, ...
                                               nRotations, ...
                                               useNormXCorr, ...
                                               corner_radius, ...
                                               excluded_xyrs) ;

% Decare victory
fprintf('Detected led indicator.\n');

end  % function

