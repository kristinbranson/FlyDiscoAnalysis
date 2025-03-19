function warn_if_protocol_is_longer_than_video(expdir, settingsdir, analysis_protocol, do_run)
% Throws an error if this is an optogenetic experiment and the protocol is
% longer than the video.  Otherwise, exits normally.  This reads all the info
% it needs from the experiment folder and the analysis-protocol folder, even
% though some of that stuff is likely already in memory when this is called.

% Get locations of parameter files
datalocparamsfile = fullfile(settingsdir,analysis_protocol,'dataloc_params.txt');
dataloc_params = ReadParams(datalocparamsfile);

% Get one thing from the indicator params
indicatorparamsfile = fullfile(settingsdir,analysis_protocol,dataloc_params.indicatorparamsfilestr);
if exist(indicatorparamsfile,'file'),
  indicator_params = loadIndicatorParams(indicatorparamsfile) ;
  isOptogeneticExp = logical(indicator_params.OptogeneticExp) ;
else
  isOptogeneticExp = false ;
end

% This check only makes sense for opto experiments
if ~isOptogeneticExp ,
  return
end

% Read in registration params
analysis_protocol_folder_path = fullfile(settingsdir, analysis_protocol) ;
registrationparamsfile = fullfile(analysis_protocol_folder_path,dataloc_params.registrationparamsfilestr);
if ~exist(registrationparamsfile,'file')
  error('Registration params file %s does not exist',registrationparamsfile);
end
raw_registration_params = ReadParams(registrationparamsfile);
registration_params = modernizeRegistrationParams(raw_registration_params) ;

% Load trajectories
ctraxfile = determine_downstream_trx_file_path(expdir, dataloc_params, do_run) ;
moviefile = fullfile(expdir,dataloc_params.moviefilestr);
[trx,~,succeeded,timestamps] = load_tracks(ctraxfile,moviefile,'annname','');
if ~succeeded,
  error('Could not load trajectories from file %s',ctraxfile);
end

% Get the frame interval for the movie, the dt, and hence the duration
[are_timestamps_reliable, fallback_dt] = getFallbackDtIfNeeded(registration_params, timestamps, trx) ;
dt = fif(are_timestamps_reliable, median(diff(timestamps)), fallback_dt) ;  % frame interval, in seconds
nframes = numel(timestamps) ;
video_duration = dt*nframes ;

% Read in the LED protocol, compute the total protocol duration
ledprotocolfile = fullfile(expdir,dataloc_params.ledprotocolfilestr);
raw_protocol = loadAnonymous(ledprotocolfile) ;
protocol = downmixProtocolIfNeeded(raw_protocol, indicator_params) ;
  % We downmix the (possibly RGB) protocol as an expedient, because we want to
  % take an RGB protocol as input.  When we enrich support for RGB protocols, we
  % should change this.
% See the comments in indicatorframes_from_protocol.m for details about how
% the protocol structure works, and what it means.
duration_from_step_index = protocol.duration/1000 ;  % seconds
active_duration_from_step_index = ...
  protocol.delayTime + protocol.iteration .* (protocol.pulseNum .* protocol.pulsePeriodSP + protocol.offTime)/1000 ;  % seconds
  % For each i, active_duration_from_step_index(i) should be <=
  % duration_from_step_index(i)
if any(active_duration_from_step_index>duration_from_step_index)
  warning('Active duration of at least one protocol step is greater than the step duration') ;
  duration_from_step_index  %#ok<NOPRT> 
  active_duration_from_step_index  %#ok<NOPRT> 
end
total_protocol_duration = sum(duration_from_step_index) ;

% Finally, compare the protocol duration to the video duration and error if
% the protocol is longer.
if total_protocol_duration > video_duration ,
  warning('According to the protocol file, the protocol is longer in duration (%g s) than the video (%g s)', ...
          total_protocol_duration, ...
          video_duration) ;
end
