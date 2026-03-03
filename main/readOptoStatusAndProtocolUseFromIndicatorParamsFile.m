function [isOptogeneticExp, doesUseProtocolDotMat, indicator_params_or_empty] = ...
  readOptoStatusAndProtocolUseFromIndicatorParamsFile(settingsdir, analysis_protocol, indicatorparamsfilestr) 
% Read a couple pivotal quantities from the indicator params file.

indicatorparamsfile = fullfile(settingsdir, analysis_protocol, indicatorparamsfilestr) ;

% Get a few things from the indicator params
if exist(indicatorparamsfile,'file'),
  indicator_params_or_empty = loadIndicatorParams(indicatorparamsfile) ;
  isOptogeneticExp = logical(indicator_params_or_empty.OptogeneticExp) ;
  doesUseProtocolDotMat = logical(indicator_params_or_empty.doesUseProtocolDotMat) ;
else
  indicator_params_or_empty = [] ;
  isOptogeneticExp = false ;
  doesUseProtocolDotMat = false ;  % shouldn't matter
end

end  % function
