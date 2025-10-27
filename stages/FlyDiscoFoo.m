function FlyDiscoFoo(expdir, varargin)
% Example of a FlyDisco stage function

% Deal with optional arguments
[analysis_protocol, settingsdir, datalocparamsfilestr, forcecompute, debug, do_run] = ...
  myparse(varargin,...
  'analysis_protocol','current_bubble',...
  'settingsdir', default_settings_folder_path(), ...
  'datalocparamsfilestr','dataloc_params.txt', ...
  'forcecompute', false, ...
  'debug',false, ...
  'do_run', []) ;

% Read in the data file locations
analysis_protocol_folder_path = fullfile(settingsdir, analysis_protocol) ;
datalocparamsfile = fullfile(analysis_protocol_folder_path, datalocparamsfilestr) ;
dataloc_params = ReadParams(datalocparamsfile) ;

% Load in the parameters file
fooParamsFilePath = fullfile(analysis_protocol_folder_path, dataloc_params.fooparamsfilestr) ;
fooParams = ReadParams(fooParamsFilePath) ;

% Determine the output file path
outputFileName = dataloc_params.foooutputfilestr ;
outputFilePath = fullfile(expdir, outputFileName) ;

% 'Touch' the output file
% touch(outputFilePath) ;  % Create an empty output file

% Compute a zany statistic from the tracks
trxFilePath = fullfile(expdir, dataloc_params.trxfilestr) ;
metatrx = load(trxFilePath) ;
trx = metatrx.trx ;
crazy_statistic_from_tracklet_index = arrayfun(@(tracklet)(mean(tan(tracklet.x))), trx, 'UniformOutput', true) ;
crazy_statistic = prod(crazy_statistic_from_tracklet_index) ;

% Write the crazy statistic to the output file
fo = file_object(outputFilePath, 'wt') ;
fo.fprintf('crazy: %.17g\n', crazy_statistic) ;

end  % function
