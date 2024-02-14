function registration_params_cell = marshallRegistrationParamsForDetectRegistrationMarks(registration_params, expdir, dataloc_params) 
% Collect things from registration_params, make a cell array suitable for
% passing to detectRegistrationMarks().
% This is a pure function.

fnsignore = intersect(fieldnames(registration_params),...
  {'minFliesLoadedTime','maxFliesLoadedTime','extraBufferFliesLoadedTime','usemediandt','doTemporalRegistration','OptogeneticExp', ...
   'LEDMarkerType','maxDistCornerFrac_LEDLabel','doTemporalTruncation','maxFlyTrackerNanInterpFrames', 'bkgdNSampleFrames'});
registration_params_cell = struct2paramscell(rmfield(registration_params,fnsignore));

% file to save image to
if isfield(dataloc_params,'registrationimagefilestr'),
  registration_params_cell(end+1:end+2) = {'imsavename',fullfile(expdir,dataloc_params.registrationimagefilestr)};
end

end
