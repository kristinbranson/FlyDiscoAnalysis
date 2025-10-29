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
% [~,outputFileBase,~] = fileparts(outputFileName) ;

% 'Touch' the output file
% touch(outputFilePath) ;  % Create an empty output file

% Plot some random thing from the registered trx
trxFilePath = fullfile(expdir, dataloc_params.trxfilestr) ;
metatrx = load(trxFilePath) ;
trx = metatrx.trx ;
x = { trx.x } ;
fig = figure('color', 'w') ;
cleaner = onCleanup(@()(close(fig))) ;
plot(x{1}, 'k') ;

% "Print" that to a .png
print(fig, outputFilePath, '-dpng', '-r300') ;

end  % function
