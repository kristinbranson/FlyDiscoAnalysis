function [isOptogeneticExp, doesUseProtocolDotMat] = ...
  readOptoStatusAndProtocolUseFromIndicatorParamsFile(settingsdir, analysis_protocol, indicatorparamsfilestr) 
% Read a couple pivotal quantities from the indicator params file.

indicatorparamsfile = fullfile(settingsdir, analysis_protocol, indicatorparamsfilestr) ;

% Get a few things from the indicator params
if exist(indicatorparamsfile,'file'),
  indicator_params = loadIndicatorParams(indicatorparamsfile) ;
  isOptogeneticExp = logical(indicator_params.OptogeneticExp) ;
  doesUseProtocolDotMat = logical(indicator_params.doesUseProtocolDotMat) ;
else
  isOptogeneticExp = false ;
  doesUseProtocolDotMat = false ;  % shouldn't matter
end

end  % function
