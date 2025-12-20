function error_if_protocol_stim_num_notequal_detected(expdir, settingsdir, analysis_protocol, varargin)
% Throws an error if this is an optogenetic experiment and the number of
% stimuli from the protocol does not match the number of stimuli detected
% from the indicator LEDs. 

% Parse the arguments
indicatordata = myparse(varargin,'indicatordata',[]);

% Get locations of parameter files
datalocparamsfile = fullfile(settingsdir,analysis_protocol,'dataloc_params.txt');
dataloc_params = ReadParams(datalocparamsfile);

% Get a few things from the indicator params
[isOptogeneticExp, doesUseProtocolDotMat] = ...
  readOptoStatusAndProtocolUseFromIndicatorParamsFile(settingsdir, analysis_protocol, dataloc_params.indicatorparamsfilestr) ;

% This check only makes sense for opto experiments with protocol
if ~(isOptogeneticExp && doesUseProtocolDotMat)
  return
end

% Read in the LED protocol, compute the total protocol duration
ledprotocolfile = fullfile(expdir,dataloc_params.ledprotocolfilestr);
if ~exist(ledprotocolfile, 'file') 
  warning('LED protocol file %s does not exist, so not comparing number of LED pulses to a protocol', ledprotocolfile) ;
  return
end
raw_protocol = loadAnonymous(ledprotocolfile) ;
protocol = downmixProtocolIfNeeded(raw_protocol, indicator_params) ;

% If indictor not passed in, load detected stim data
if isempty(indicatordata)
  indicatordatafile = fullfile(expdir,dataloc_params.indicatordatafilestr);
  load(indicatordatafile,'indicatorLED');
else
  indicatorLED = indicatordata.indicatorLED;
end
downmixed_indicatorLED = downmix_indicatorLED(indicatorLED) ;

% Perform the test, error if it fails
nprotocolstim = distinct_pulse_count_from_single_channel_protocol(protocol, indicator_params.pad) ;
ndetectedstim = numel(downmixed_indicatorLED.startframe) ;
if nprotocolstim ~= ndetectedstim
  error('The number of detected stimuli (%g) does not equal the expected number of stimuli (%g) from the protocol file', ...
        ndetectedstim, ...
        nprotocolstim) ;
end
