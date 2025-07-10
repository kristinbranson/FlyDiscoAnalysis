function [stride, lastFrameIndex] = determineStrideAndLastFrameIndexForLedMaxImage(protocol, dt, nframes)
% To make the "mex LED image", a max will be taken of a bunch of movie frames.
% These sampled frames start with the first video frame, and will be sampled
% every stride'th frame, up to frame frameToLedPulse.  This function computes
% good values for stride and frameToLedPulse.

% Compute rough values for frameToLedPulse and stride
if ~isempty(protocol)
  % Determine which step is the first with non-zero intensity
  isActiveFromStepIndex = (protocol.intensity~=0) ;  % Steps can have intensity of zero.  We don't want those.
  firstActiveStepIndex = find(isActiveFromStepIndex,1);
  if isempty(firstActiveStepIndex)
    error('No steps with LED intensity greater than zero')
  end

  % Compute the time of the end of the first active "step" in the protocol.
  durationBeforeFirstActiveStep = sum(protocol.duration(1:firstActiveStepIndex-1))/1000 ;  % seconds
  delay = protocol.delayTime(firstActiveStepIndex) ;  % seconds
  pulsePeriod = protocol.pulsePeriodSP(firstActiveStepIndex)/1000 ;  % seconds
  pulseCount = protocol.pulseNum(firstActiveStepIndex) ;  % pure
  timeOfEndOfFirstActivePulseTrain = durationBeforeFirstActiveStep + delay + pulsePeriod*pulseCount ;  % seconds

  % sets end range in which to find LED on
  movieDuration = dt*nframes ;  % seconds
  if timeOfEndOfFirstActivePulseTrain > movieDuration ,
    error('According to the protocol, the first active pulse train ends at t = %g s.  This is after the end of the video (duration = %g s).', ...
          timeOfEndOfFirstActivePulseTrain, ...
          movieDuration) ;
  end
  timeOfEndOfFirstActivePulseTrainInFrameIntervals = timeOfEndOfFirstActivePulseTrain/dt ;
  roughStride = max(1, round(pulsePeriod*pulseCount/(2*dt))) ;
else
  % reasonable guess if no protocol
  timeOfEndOfFirstActivePulseTrainInFrameIntervals = nframes/3 ;
  roughStride = 10 ;
end

% Put the final touches on things
stride = min(roughStride, round(timeOfEndOfFirstActivePulseTrainInFrameIntervals/100)) ;
lastFrameIndex = max(round(timeOfEndOfFirstActivePulseTrainInFrameIntervals), min(200, nframes)) ;
