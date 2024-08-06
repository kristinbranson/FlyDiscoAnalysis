function indicatorframes = indicatorframes_from_protocol(protocol)
% Given the protocol, determine which stimuli to use for each video snippet.
% indicatorframes(snippet_index) gives the stimulus index for that snippet.
% Here a "stimulus" is a single "pulse train" as described below.
% "indicatorframes" could maybe also be named train_index_from_snippet_index,
% if one was being long-winded, I think  -- ALT, 2024-03-22
%
% The protocol describes a pattern of pulses to be delivered.  These pulses
% are organized into "pulse trains", each train an evenly-spaced series of
% individual pulses.  A single "step" of the protocol specifies a
% "supertrain": a series of pulse trains.
%
% An example protocol structure looks like this:
%
% protocol =
%
%   struct with fields:
%
%           duration: [22×1 double] (ms)
%          intensity: [22×1 double]
%       pulseWidthSP: [22×1 double] (ms)
%      pulsePeriodSP: [22×1 double] (ms)
%           pulseNum: [22×1 double] (pure, a count)
%            offTime: [22×1 double] (ms)
%          delayTime: [22×1 double] (s)
%          iteration: [22×1 double] (pure, a count)
%
% This protocol contains 22 steps, i.e. 22 supertrains.  Step i is of duration
% duration(i).  Each step begins with a delay.  The delay for step i is
% delayTime(i).  This is followed by a pulse train.  The number of pulses in
% step i is given by pulseNum(i).  The duration of each pulse in step i is
% given by pulseWidthSP(i) (in ms).  The time between pulse starts in step i
% is given by pulsePeriodSP(i) (in ms).  At the end of each pulse train is an
% "off time", where no pulses occur.  The duration of each off time in step i
% is given by offTime(i).  The pulse train consists of the evenly-spaced
% pulses plus the off time at the end of each pulse train.  Thus each pulse
% train in step i is of duration pulseNum(i)*pulsePeriodSP(i)/1000+offTime(i).
% Step i contains iteration(i) pulse trains.  Thus step i consists of a delay
% of duration delayTime(i), then a series of iteration(i) pulse trains.  The
% "active" duration of step i is thus:
%
%   delayTime(i) + iteration(i) * (pulseNum(i)*pulsePeriodSP(i)+offTime(i))/1000
%
% For each step, the active duration should be less than or equal to the
% duration.  Any leftover time at the end of each step will contain no pulses.
% The total duration of the protocol is the sum of the step durations.
%
% The amplitude of all the pulses in step i is intensity(i).

step_count = numel(protocol.stepNum);

% Figure out some things about the mapping from train indices to step indices
preceding_train_count_from_step_index = [0 ; cumsum(protocol.iteration(1:end-1))] ;  % The number of trains that come before each step
first_train_index_from_step_index = preceding_train_count_from_step_index + 1 ;  % The train index of the first train in each step

% Figure out which steps, if any, are "blank" (i.e. intensity==0)
intensity_from_step_index = protocol.intensity ;
is_nonblank_from_step_index = ( intensity_from_step_index>0 ) ;
nonblank_step_count = sum(is_nonblank_from_step_index) ;
step_index_from_nonblank_step_index = find(is_nonblank_from_step_index) ;

% How snippets are chosen depends on the step count
if nonblank_step_count == 0 ,
  error('There are zero non-blank steps') ;
elseif nonblank_step_count == 1,
  step_index = step_index_from_nonblank_step_index(1) ;
  first_train_index = first_train_index_from_step_index(step_index) ;  % index of first train in the step
  train_count = protocol.iteration(step_index) ;  % Number of pulse trains in the step
  train_stride = ceil(train_count/6);  % take every stride'th iteration of the stimulus
  indicatorframes = (first_train_index-1) + 1:train_stride:train_count ;
  step_index_from_snippet_index = ones(size(indicatorframes)) ;
elseif nonblank_step_count <= 3,
  % In this case we do two snippets for each nonblank step, one for the first
  % train of the step, one for the last train of the step.
  % We assume all steps have more 2 or more iterations/trains
  if ~all(protocol.iteration>=2)
    error('There are two or three nonblank steps, and at least one step with iterations < 2. User needs to make a specific ctrax results movie param file')
  end
  snippet_count = 2*nonblank_step_count ;
  indicatorframes = zeros(1,snippet_count);
  step_index_from_snippet_index = zeros(1,snippet_count) ;
  for nonblank_step_index = 1:nonblank_step_count ,
    step_index = step_index_from_nonblank_step_index(nonblank_step_index) ;
    first_train_index = first_train_index_from_step_index(step_index) ;  % Index of first train in the step
    train_count = protocol.iteration(step_index);
    first_snippet_index = 2*(nonblank_step_index-1)+1 ;
    second_snippet_index = 2*(nonblank_step_index-1)+2 ;
    indicatorframes(first_snippet_index) = first_train_index ;
    indicatorframes(second_snippet_index) = first_train_index+train_count-1 ;
    step_index_from_snippet_index(first_snippet_index:second_snippet_index) = step_index ;
  end
else
  % If there are many steps, do one snippet for the first train in each step.
  nonblank_step_stride = ceil(step_count/6) ;
  nonblank_step_index_from_snippet_index = 1:nonblank_step_stride:nonblank_step_count ;
  step_index_from_snippet_index = step_index_from_nonblank_step_index(nonblank_step_index_from_snippet_index) ;
  raw_indicatorframes = first_train_index_from_step_index(step_index_from_snippet_index) ;  % indicatorframes, but as a col vector
  indicatorframes = raw_indicatorframes(:)' ;  % Make a row vector
end

% Double-check that none of the steps from which snippets are drawn have zero intensity
intensity_from_snippet_index = protocol.intensity(step_index_from_snippet_index) ;
% If any remaining snippets are from steps with zero intensity, raise an error
if any(intensity_from_snippet_index == 0) ,
  error('Internal error: At least one of the snippets is from a step with zero intensity')
end
