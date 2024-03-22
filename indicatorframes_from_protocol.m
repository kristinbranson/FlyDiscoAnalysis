function indicatorframes = indicatorframes_from_protocol(protocol)
% Given the protocol, determine which stimuli to use for each video snippet.
% indicatorframes(snippet_index) gives the stimulus index for that snippet.
% Here a "stimulus" is a single "pulse train" as described below.
% "indicatorframes" could maybe also be named train_index_from_snippet_index,
% if one was being long-winded, I think  -- ALT, 2024-03-22

% The protocol describes a pattern of pulses to be delivered.  These pulses
% are organized into "pulse trains", each train an evenly-spaced series of
% individual pulses.  A single "step" of the protocol specifies a
% "supertrain": a series of pulse trains.

% An example protocol structure looks like this:
%
% protocol =
%
%   struct with fields:
%
%          intensity: [22×1 double]
%       pulseWidthSP: [22×1 double]
%      pulsePeriodSP: [22×1 double]
%           pulseNum: [22×1 double]
%            offTime: [22×1 double]
%          delayTime: [22×1 double]
%          iteration: [22×1 double]
%
% This protocol contains 22 steps, i.e. 22 supertrains.  Each step
% begins with a delay.  The delay for step i is delayTime(i).  This is
% followed by a pulse train.  The number of pulses in step i is given by
% pulseNum(i).  The duration of each pulse in step i is given by
% pulseWidthSP(i) (in ms).  The time between pulse starts in step i is given
% by pulsePeriodSP(i) (in ms).  At the end of each pulse train is an "off time", where
% no pulses occur.  The duration of each off time in step i is given by
% offTime(i).  The pulse train consists of the evenly-spaced pulses plus the
% off time at the end of each pulse train.  Thus each pulse train in step i is
% of duration pulseNum(i)*pulsePeriodSP(i)/1000+offTime(i).  Step i contains
% iteration(i) pulse trains.  Thus step i consists of a delay of duration
% delayTime(i), then a series of iteration(i) pulse trains.  The duration of
% step i is thus:
%
%   delayTime(i) + iteration(i) * (pulseNum(i)*pulsePeriodSP(i)/1000+offTime(i))
%
% The total duration of the protocol is the sum of the step durations.
%
% The amplitude of all the pulses in step i is intensity(i).

step_count = numel(protocol.stepNum);
% if there's only 1 step
if step_count == 1,
  train_count = protocol.iteration;
  train_stride = ceil(train_count/6);     % take every stride'th iteration of the stimulus
  indicatorframes = 1:train_stride:train_count;
  step_index_from_snippet_index = ones(size(indicatorframes)) ;
elseif step_count <= 3,
  % assumes steps have more 2 or more iterations/trains
  if ~all([protocol.iteration] > 1)
    error('Step with iteration < 2. User needs to make a specific ctrax results movie param file')
  end
  snippet_count = 2*step_count ;
  indicatorframes = zeros(1,snippet_count);
  step_index_from_snippet_index = zeros(1,snippet_count) ;
  for step_index = 1:step_count,
    train_count = protocol.iteration(step_index);
    if step_index==1,
      indicatorframes(1)=1;
      indicatorframes(2)=train_count;
      step_index_from_snippet_index(1:2) = 1 ;
    elseif step_index == 2
      indicatorframes(3) = indicatorframes(2)+1;
      indicatorframes(4) = indicatorframes(2)+train_count;
      step_index_from_snippet_index(3:4) = 2 ;
    elseif step_index == 3
      indicatorframes(5) = indicatorframes(4)+1;
      indicatorframes(6) = indicatorframes(4)+train_count;
      step_index_from_snippet_index(5:6) = 3 ;
    end
  end
else
  % If there are many steps, do one snippet for the first train in each step.
  step_stride = ceil(step_count/6) ;
  step_index_from_snippet_index = 1:step_stride:step_count ;
  preceding_train_count_from_step_index = [0 ; cumsum(protocol.iteration(1:end-1))] ;  % The number of trains that come before each step
  first_train_index_from_step_index = preceding_train_count_from_step_index + 1 ;  % The train index of the first train in each step
  raw_indicatorframes = first_train_index_from_step_index(step_index_from_snippet_index) ;  % indicatorframes, but as a col vector
  indicatorframes = raw_indicatorframes(:)' ;  % Make a row vector
end

% Make sure none of the steps from which snippets are drawn have zero intensity
snippet_count = numel(indicatorframes) ;
intensity_from_snippet_index = protocol.intensity(step_index_from_snippet_index) ;
% If the last snippet is from a step with zero intensity, just drop that
% snippet
if intensity_from_snippet_index(snippet_count) == 0
  indicatorframes(snippet_count) = [] ;
  intensity_from_snippet_index(snippet_count) = [] ;
end
% If any remaining snippets are from steps with zero intensity, raise an error
if any(intensity_from_snippet_index == 0) ,
  error('Step with intensity = 0 in middle of experiment. User needs to make specific ctrax results movie params file')
end
