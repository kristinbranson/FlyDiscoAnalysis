function [ledMaxImage, ledTemplate] = ...
  determineLedMaxImageAndLedTemplate(expdir, ...
                                     analysis_protocol_folder_path, ...
                                     dataloc_params, ...
                                     indicator_params, ...
                                     registration_data, ...
                                     protocol, ...
                                     dt, ...
                                     rigId)

% Extract an image from the movie that will be good for finding the LED.  Also
% read the LED template image from the analysis-protocol folder.

% protocol should be in ChR (single LED) format (i.e. should have .intensity
% field, not .Rintensity, etc.)

if ~isfield(indicator_params,'LEDMarkerType') ,
  error('No LEDMarkerType field in indicator params') ;
end

% Get the movie header and a frame-reader function
moviefilename = fullfile(expdir,dataloc_params.moviefilestr);
[readfcn,~,~,header] = get_readframe_fcn(moviefilename);
nframes = header.nframes ;

% Determine the stride for sampling movie frames to make the ledMaxImage
[stride, lastFrameIndex] = determineStrideAndLastFrameIndexForLedMaxImage(protocol, dt, nframes) ;
firstPassLedMaxImage = sampleFramesForMaximumImage(readfcn, 1, stride, lastFrameIndex) ;

% Determine the path to the LED template image
template_file_path = determine_template_or_mask_file_path(indicator_params.LEDMarkerType, rigId, analysis_protocol_folder_path) ;

%
% Check if something of similar brightness to the LED template is present in
% the ledMaxImage.  If not, make a new ledMaxImage.
%
ledTemplate = double(imread(template_file_path)) ;
binWidth = 25;
binEdges = (0:binWidth:255) ;  % [ 0    25    50    75   100   125   150   175   200   225   250 ], length==11
templatePixelCountFromBinIndex = histcounts(ledTemplate,binEdges);  % length==10, note does not count pixel values >250
brightPixelCountInTemplate = sum(templatePixelCountFromBinIndex(end-1:end));  % Total number of pixels in the top two bins
firstFrame = readfcn(1) ;
diffImage = firstPassLedMaxImage - firstFrame ;

[nr,nc,~] = size(firstPassLedMaxImage) ;
%[xgrid, ygrid] = meshgrid(1:nc, 1:nr);
%outsideArenaMask = ((xgrid-registration_data.circleCenterX).^2 + (ygrid-registration_data.circleCenterY).^2) >= registration_data.circleRadius.^2;
outsideArenaMask = outsideArenaMaskFromRegistrationData(registration_data, nc, nr) ;
outsideArenaPixels = diffImage(outsideArenaMask);  % get the pixel values outside the arena
outsideArenaPixelCountFromBinIndex = histcounts(outsideArenaPixels,binEdges);  % make a histogram of those values
brightPixelCountOutsideArena = sum(outsideArenaPixelCountFromBinIndex(end-1:end));  % get the total number of pixels in the top two bins
if brightPixelCountOutsideArena <= brightPixelCountInTemplate
  % The number of bright pixels outside the arena in the diffimage seems to be
  % less than or equal to the number of bright pixels in the template image.
  % This suggests the LED template is not present in the maxLedImage.
  % So we make a new ledMaxImage, starting from the first movie frame, and hope
  % for better results.
  ledMaxImage = sampleFramesForMaximumImage(readfcn, 1, stride, nframes) ;
else
  ledMaxImage = firstPassLedMaxImage ;
end
