function oneStep = ProtocolExperimentSteps(protocol)

oneStep = struct();

%calculate the four point value for each step
stepCount = length(protocol.stepNum) ;
for stepIndex = 1:stepCount
  oneStep(stepIndex).NumStep = protocol.stepNum(stepIndex);
  oneStep(stepIndex).Duration = protocol.duration(stepIndex)/1000; % specify in ms but firmware expects seconds;
  oneStep(stepIndex).DelayTime = protocol.delayTime(stepIndex);
  %red light
  if isfield(protocol,'Rintensity') 
    oneStep(stepIndex).RedIntensity = protocol.Rintensity(stepIndex);
    oneStep(stepIndex).RedPulseWidth = protocol.RpulseWidth(stepIndex);
    oneStep(stepIndex).RedPulsePeriod = protocol.RpulsePeriod(stepIndex);
    oneStep(stepIndex).RedPulseNum = protocol.RpulseNum(stepIndex);
    oneStep(stepIndex).RedOffTime = protocol.RoffTime(stepIndex);
    oneStep(stepIndex).RedIteration = protocol.Riteration(stepIndex);
  else
    oneStep(stepIndex).RedIntensity = 0 ;
    oneStep(stepIndex).RedPulseWidth = 0 ;
    oneStep(stepIndex).RedPulsePeriod = 0 ;
    oneStep(stepIndex).RedPulseNum = 0 ;
    oneStep(stepIndex).RedOffTime = 0 ;
    oneStep(stepIndex).RedIteration = 0 ;
  end    
  %green light
  if isfield(protocol,'Gintensity') 
    oneStep(stepIndex).GrnIntensity = protocol.Gintensity(stepIndex);
    oneStep(stepIndex).GrnPulseWidth = protocol.GpulseWidth(stepIndex);
    oneStep(stepIndex).GrnPulsePeriod = protocol.GpulsePeriod(stepIndex);
    oneStep(stepIndex).GrnPulseNum = protocol.GpulseNum(stepIndex);
    oneStep(stepIndex).GrnOffTime = protocol.GoffTime(stepIndex);
    oneStep(stepIndex).GrnIteration = protocol.Giteration(stepIndex);
  else
    oneStep(stepIndex).GrnIntensity = 0 ;
    oneStep(stepIndex).GrnPulseWidth = 0 ;
    oneStep(stepIndex).GrnPulsePeriod = 0 ;
    oneStep(stepIndex).GrnPulseNum = 0 ;
    oneStep(stepIndex).GrnOffTime = 0 ;
    oneStep(stepIndex).GrnIteration = 0 ;
  end    
  %blue light  
  if isfield(protocol,'Bintensity') 
    oneStep(stepIndex).BluIntensity = protocol.Bintensity(stepIndex);
    oneStep(stepIndex).BluPulseWidth = protocol.BpulseWidth(stepIndex);
    oneStep(stepIndex).BluPulsePeriod = protocol.BpulsePeriod(stepIndex);
    oneStep(stepIndex).BluPulseNum = protocol.BpulseNum(stepIndex);
    oneStep(stepIndex).BluOffTime = protocol.BoffTime(stepIndex);
    oneStep(stepIndex).BluIteration = protocol.Biteration(stepIndex);
  else
    oneStep(stepIndex).BluIntensity = 0 ;
    oneStep(stepIndex).BluPulseWidth = 0 ;
    oneStep(stepIndex).BluPulsePeriod = 0 ;
    oneStep(stepIndex).BluPulseNum = 0 ;
    oneStep(stepIndex).BluOffTime = 0 ;
    oneStep(stepIndex).BluIteration = 0 ;
  end
end

