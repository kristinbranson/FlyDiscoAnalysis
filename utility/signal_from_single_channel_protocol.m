function [y, t] = signal_from_single_channel_protocol(protocol, dt)
% Given a single-channel protocol, generate the corresponding LED intensity vs
% time signel.
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

% Deal with arguments
if ~exist('dt', 'var') || isempty(dt) ,
  dt = 0.001 ;  % s
end

% Unpack the protocol into individual variables, converting all times to
% seconds
duration_from_step_index = protocol.duration/1000 ;  % ms->s
intensity_from_step_index = protocol.intensity ;  % arbs
pulse_duration_from_step_index = protocol.pulseWidthSP/1000 ;  % ms->s
pulse_period_from_step_index = protocol.pulsePeriodSP/1000 ;  % ms->s
pulse_count_from_step_index = protocol.pulseNum ;  % pure, a count
inter_train_interval_from_step_index = protocol.offTime/1000 ;  % ms->s
train_delay_from_step_index = protocol.delayTime ;  % s
train_count_from_step_index = protocol.iteration ;  % pure, a count

% Get the number of steps, aka number of supertrains
step_count = numel(duration_from_step_index) ;

% Compute some derived quantities we'll need
offset_time_from_step_index = cumsum(duration_from_step_index) ;
onset_time_from_step_index = vertcat(0, offset_time_from_step_index(1:end-1)) ;
%pulse_dc_from_step_index = pulse_duration_from_step_index ./ pulse_period_from_step_index ;
train_duration_from_step_index = pulse_period_from_step_index .* pulse_count_from_step_index ;
train_period_from_step_index = train_duration_from_step_index + inter_train_interval_from_step_index ;
%train_dc_from_step_index = train_duration_from_step_index ./ train_period_from_step_index ;
% Each step has a "window" that starts after the train delay, and is the
% duration of the train count times the train period.
window_duration_from_step_index = train_period_from_step_index .* train_count_from_step_index ;
window_onset_from_step_index = onset_time_from_step_index + train_delay_from_step_index ;

% Create the signal timeline
total_duration = sum(duration_from_step_index) ;  % s
n = round(total_duration/dt) ;  % number of elements in t, y
t = dt/2+dt*(0:(n-1))' ;  % 1xn, seconds, offset by half a sample to hopefully avoid samples at signal edges

% For each step, generate a signal for that step.  Final signal is the sum of
% these.
y = zeros(size(t)) ;
for step_index = 1 : step_count ,
  % Extract needed params for this step
  window_onset_time = window_onset_from_step_index(step_index) ;
  window_duration = window_duration_from_step_index(step_index) ;
  pulse_period = pulse_period_from_step_index(step_index) ;
  pulse_duration = pulse_duration_from_step_index(step_index) ;
  train_period = train_period_from_step_index(step_index) ;
  train_duration = train_duration_from_step_index(step_index) ;
  intensity = intensity_from_step_index(step_index) ;  

  % Compute the signal just for this step
  t_window = t - window_onset_time ;  % time relative to onset of current step's window
  y_this_step = window(t_window, window_duration) .* ...
                square_wave(t_window, train_period, train_duration) .* ...
                square_wave(t_window, pulse_period, pulse_duration) .* ...
                intensity ;

  % Add this step's signal to the running total
  y = y + y_this_step ;
end

end  % function



function y = square_wave(t, T, dur)
% A square wave with amplitude 1, period T, and duty-cycle dur/T.
% Sampled at times in t.  The first rising edge will be at t==0.  A square
% wave is also known as a pulse train.
phase = mod(t,T) ;
y = double(phase < dur) ;
end  % function



function y = window(t, T)
% A "window" function that is one on 0<=t<T, otherwise zero.
y = double( (0 <= t) .* (t < T ) ) ;
end  % function
