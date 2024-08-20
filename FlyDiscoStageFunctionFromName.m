function result = FlyDiscoStageFunctionFromName()

% Return a scalar struct mapping FlyDisco stage names to stage function handles.
% The order is important---this is the order they will be run in the pipeline.

result = struct() ;
result.automaticchecksincoming = @FlyDiscoAutomaticChecksIncoming ;
result.flytracking = @FlyTrackerWrapperForFlyDisco ;
result.addpflies = @FlyDiscoAddPFlies ;
result.registration = @FlyDiscoRegisterTrx ;
result.ledonoffdetection = @FlyDiscoDetectIndicatorLedOnOff ;
result.sexclassification = @FlyDiscoClassifySex ;
result.apt = @FlyDiscoAPTTrackWrapper ;
result.computeperframefeatures = @FlyDiscoComputePerFrameFeatures ;
result.computehoghofperframefeatures = @FlyDiscoComputeHOGHOFPerFrameFeatures ;
result.jaabadetect = @FlyDiscoJAABADetect ;
result.locomotionmetrics = @FlyDiscoComputeLocomotionMetrics ;
result.computeperframestats = @FlyDiscoComputePerFrameStats ;
result.plotperframestats = @FlyDiscoPlotPerFrameStats ;
result.makectraxresultsmovie = @FlyDiscoMakeCtraxResultsMovie ;
result.makeaptresultsmovie = @FlyDiscoMakeAPTResultsMovie ;
result.automaticcheckscomplete = @FlyDiscoAutomaticChecksComplete ;

end
