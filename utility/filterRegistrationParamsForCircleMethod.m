function registration_params_cell = filterRegistrationParamsForCircleMethod(registration_params) 
% Collect things from registration_params, make a cell array suitable for
% passing to determineRegistrationTransformForCircleMethodCore().
% This is a pure function.

fnsignore = intersect(fieldnames(registration_params),...
  {'minFliesLoadedTime','maxFliesLoadedTime','extraBufferFliesLoadedTime','usemediandt','doTemporalRegistration','OptogeneticExp', ...
   'LEDMarkerType','maxDistCornerFrac_LEDLabel','doTemporalTruncation','maxFlyTrackerNanInterpFrames', 'bkgdNSampleFrames','featureRadius',...
   'bowlMarkerType','maxDistCornerFrac_BowlLabel','method'});
registration_params_cell = struct2paramscell(rmfield(registration_params,fnsignore));

end
