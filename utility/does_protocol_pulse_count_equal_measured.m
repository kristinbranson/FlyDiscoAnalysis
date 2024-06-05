function [do_stim_counts_match, message] = does_protocol_pulse_count_equal_measured(expdir, settingsdir, analysis_protocol, varargin)
% Throws an error if this is an optogenetic experiment and the number of
% stimuli from the protocol does not match the number of stimuli detected
% from the indicator LEDs.

% Parse the arguments
[indicatordata, dothrowerror] = myparse(varargin, 'indicatordata', [], 'dothrowerror', false) ;

% Get locations of parameter files
datalocparamsfile = fullfile(settingsdir,analysis_protocol,'dataloc_params.txt');
dataloc_params = ReadParams(datalocparamsfile);

% Don't need this check---this is called from a location that is only reached
% if it's an optogenetic experiment.
%
% % Get one thing from the indicator params
% indicatorparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.indicatorparamsfilestr);
% if exist(indicatorparamsfile,'file'),
%   raw_indicator_params = ReadParams(indicatorparamsfile);
%   indicator_params = modernizeIndicatorParams(raw_indicator_params) ;
%   isOptogeneticExp = logical(indicator_params.OptogeneticExp) ;
% else
%   isOptogeneticExp = false ;
% end
%
% % This check only makes sense for opto experiments
% if ~isOptogeneticExp ,
%   return
% end

% read in protocol file if it exists
ledprotocolfile = fullfile(expdir,dataloc_params.ledprotocolfilestr);
indicatordatafile = fullfile(expdir,dataloc_params.indicatordatafilestr);

if exist(ledprotocolfile,'file') && exist(indicatordatafile,'file')
    raw_protocol = loadAnonymous(ledprotocolfile) ;
    protocol = downmixProtocolIfNeeded(raw_protocol) ;


    % If indictor not passed in, load detected stim data
    if isempty(indicatordata)
        indicatordatafile = fullfile(expdir,dataloc_params.indicatordatafilestr);
        load(indicatordatafile,'indicatorLED');
    else
        indicatorLED = indicatordata.indicatorLED;
    end

    % Perform the test, error if it fails
    nprotocolstim = distinct_pulse_count_from_single_channel_protocol(protocol) ;
    ndetectedstim = numel(indicatorLED.startframe) ;
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
        end
    end
else
    do_stim_counts_match = false;
    message = sprintf('Missing LED files: could not check stimuli number');
    if dothrowerror ,
        error(message) ;  %#ok<SPERR>
    end
end