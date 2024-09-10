function [do_stim_counts_match, message] = does_protocol_pulse_count_equal_measured(expdir, settingsdir, analysis_protocol, varargin)
% Throws an error if the number of stimuli from the protocol does not match
% the number of stimuli detected from the indicator LEDs.  Should only be
% called for an optogenetic experiment.

% Parse the arguments
[indicatordata, dothrowerror] = myparse(varargin, 'indicatordata', [], 'dothrowerror', false) ;

% Get locations of parameter files
datalocparamsfile = fullfile(settingsdir,analysis_protocol, 'dataloc_params.txt') ;
dataloc_params = ReadParams(datalocparamsfile) ;

% Load the indicator params
indicatorparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.indicatorparamsfilestr);
if exist(indicatorparamsfile,'file'),
  raw_indicator_params = ReadParams(indicatorparamsfile);
  indicator_params = modernizeIndicatorParams(raw_indicator_params) ;
  isOptogeneticExp = logical(indicator_params.OptogeneticExp) ;
else
  isOptogeneticExp = false ;
end

% This check only makes sense for opto experiments
if ~isOptogeneticExp ,
  return
end

% Check that LED protocol file exists
ledprotocolfile = fullfile(expdir,dataloc_params.ledprotocolfilestr);
if ~exist(ledprotocolfile,'file') ,
    do_stim_counts_match = false;
    message = sprintf('Could not check number of stimuli: Missing LED protocol file %s', ledprotocolfile) ;
    if dothrowerror ,
        error(message) ;  %#ok<SPERR>
    else
        return
    end
end

% Load the LED protocol file, downmix it to a single-channel protocol
raw_protocol = loadAnonymous(ledprotocolfile) ;
protocol = downmixProtocolIfNeeded(raw_protocol, indicator_params) ;

% If indicator not passed in, load indicator data (as gleaned from the movie)
if isempty(indicatordata)
    % Make sure the indicator data file exists
    indicatordatafile = fullfile(expdir,dataloc_params.indicatordatafilestr) ;
    if ~exist(indicatordatafile,'file') ,
        do_stim_counts_match = false;
        message = sprintf('Could not check number of stimuli: Missing LED indicator data file %s', indicatordatafile) ;
        if dothrowerror ,
            error(message) ;  %#ok<SPERR>
        else
            return
        end
    end
    % Load the indicator data file
    load(indicatordatafile,'indicatorLED');
else
    indicatorLED = indicatordata.indicatorLED;
end
downmixed_indicatorLED = downmix_indicatorLED(indicatorLED) ;

% Perform the test, error if it fails
nprotocolstim = distinct_pulse_count_from_single_channel_protocol(protocol, indicator_params.pad) ;
ndetectedstim = numel(downmixed_indicatorLED.startframe) ;
do_stim_counts_match = ( nprotocolstim == ndetectedstim ) ;
if do_stim_counts_match ,
    message = sprintf('The number of detected stimuli (%g) is equal to the expected number of stimuli (%g) from the protocol file', ...
                      ndetectedstim, ...
                      nprotocolstim) ;
else
    message = sprintf('The number of detected stimuli (%g) does not equal the expected number of stimuli (%g) from the protocol file', ...
                      ndetectedstim, ...
                      nprotocolstim) ;
    if dothrowerror ,
        error(message) ;  %#ok<SPERR>
    else
        return
    end
end
