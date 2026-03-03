function registration_params_cell = filterRegistrationParamsForCircleMethod(registration_params) 
% Collect things from registration_params, make a cell array suitable for
% passing to determineRegistrationTransformForCircleMethodCore().
% This is a pure function.

fieldNamesToRemove = ...
  {'minFliesLoadedTime','maxFliesLoadedTime','extraBufferFliesLoadedTime','usemediandt','doTemporalRegistration','OptogeneticExp', ...
   'LEDMarkerType','maxDistCornerFrac_LEDLabel','doTemporalTruncation','maxFlyTrackerNanInterpFrames', 'bkgdNSampleFrames','featureRadius',...
   'bowlMarkerType','maxDistCornerFrac_BowlLabel','method','doesYAxisPointUp','nBowlMarkers'} ;
% rmfield() will error if you try to remove a nonexistant field, so need to
% intersect with the existing field names.
originalFieldNames = fieldnames(registration_params) ;
originalFieldNamesToRemove = intersect(originalFieldNames, fieldNamesToRemove) ;
filtered_registration_params = rmfield(registration_params,originalFieldNamesToRemove) ;
registration_params_cell = struct2paramscell(filtered_registration_params);

end
