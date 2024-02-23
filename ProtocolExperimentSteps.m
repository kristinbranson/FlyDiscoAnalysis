function oneStep = ProtocolExperimentSteps(protocol)

oneStep = struct();

%calculate the four point value for each step
for stepIndex = 1:length(protocol.stepNum)
  oneStep(stepIndex).NumStep = protocol.stepNum(stepIndex);
  oneStep(stepIndex).Duration = protocol.duration(stepIndex)/1000; % specify in ms but firmware expects seconds;
  oneStep(stepIndex).DelayTime = protocol.delayTime(stepIndex);
  %red light
  oneStep(stepIndex).RedIntensity = protocol.Rintensity(stepIndex);
  oneStep(stepIndex).RedPulseWidth = protocol.RpulseWidth(stepIndex);
  oneStep(stepIndex).RedPulsePeriod = protocol.RpulsePeriod(stepIndex);
  oneStep(stepIndex).RedPulseNum = protocol.RpulseNum(stepIndex);
  oneStep(stepIndex).RedOffTime = protocol.RoffTime(stepIndex);
  oneStep(stepIndex).RedIteration = protocol.Riteration(stepIndex);
  %green light
  oneStep(stepIndex).GrnIntensity = protocol.Gintensity(stepIndex);
  oneStep(stepIndex).GrnPulseWidth = protocol.GpulseWidth(stepIndex);
  oneStep(stepIndex).GrnPulsePeriod = protocol.GpulsePeriod(stepIndex);
  oneStep(stepIndex).GrnPulseNum = protocol.GpulseNum(stepIndex);
  oneStep(stepIndex).GrnOffTime = protocol.GoffTime(stepIndex);
  oneStep(stepIndex).GrnIteration = protocol.Giteration(stepIndex);
  %blue light
  oneStep(stepIndex).BluIntensity = protocol.Bintensity(stepIndex);
  oneStep(stepIndex).BluPulseWidth = protocol.BpulseWidth(stepIndex);
  oneStep(stepIndex).BluPulsePeriod = protocol.BpulsePeriod(stepIndex);
  oneStep(stepIndex).BluPulseNum = protocol.BpulseNum(stepIndex);
  oneStep(stepIndex).BluOffTime = protocol.BoffTime(stepIndex);
  oneStep(stepIndex).BluIteration = protocol.Biteration(stepIndex);
end
