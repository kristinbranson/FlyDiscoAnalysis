function [result, t, oneStep] = ComputeRgbStimulusProtocol(protocol)

oneStep = ProtocolExperimentSteps(protocol);

totalDuration = sum(protocol.duration)/1000; % was in ms

t = (0:0.001:totalDuration)';  % seconds
Yr = zeros(totalDuration*1000+1,1);
Yg = zeros(totalDuration*1000+1,1);
Yb = zeros(totalDuration*1000+1,1);

stepStartPnt = 1;

%calculate the four point value for each step
for stepIndex = 1:length(protocol.stepNum)   

    powerG = oneStep(stepIndex).GrnIntensity/100;
    powerR = oneStep(stepIndex).RedIntensity/100;
    powerB = oneStep(stepIndex).BluIntensity/100;
    
    LEDOnStartPnt = oneStep(stepIndex).DelayTime*1000 + stepStartPnt;
    RedOnStartPnt = LEDOnStartPnt;
    GrnOnStartPnt = LEDOnStartPnt;
    BluOnStartPnt = LEDOnStartPnt;
    
    if oneStep(stepIndex).RedIntensity > 0       
        for index = 1:oneStep(stepIndex).RedIteration
            numPntOn = oneStep(stepIndex).RedPulsePeriod*oneStep(stepIndex).RedPulseNum;
            Yr(RedOnStartPnt:RedOnStartPnt+numPntOn-1) = ones(numPntOn,1).*powerR;
            RedOnStartPnt = RedOnStartPnt + numPntOn + oneStep(stepIndex).RedOffTime-1;
        end
    end
    
    if oneStep(stepIndex).GrnIntensity > 0
        for index = 1:oneStep(stepIndex).GrnIteration
            numPntOn = oneStep(stepIndex).GrnPulsePeriod*oneStep(stepIndex).GrnPulseNum;
            Yg(GrnOnStartPnt:GrnOnStartPnt+numPntOn-1) = ones(numPntOn,1).*powerG;
            GrnOnStartPnt = GrnOnStartPnt + numPntOn + oneStep(stepIndex).GrnOffTime-1;
        end
    end
    
    if oneStep(stepIndex).BluIntensity > 0
        for index = 1:oneStep(stepIndex).BluIteration
            numPntOn = oneStep(stepIndex).BluPulsePeriod*oneStep(stepIndex).BluPulseNum;
            Yb(BluOnStartPnt:BluOnStartPnt+numPntOn-1) = ones(numPntOn,1).*powerB;
            BluOnStartPnt = BluOnStartPnt + numPntOn + oneStep(stepIndex).BluOffTime-1;
        end
    end

    stepStartPnt = stepStartPnt + oneStep(stepIndex).Duration*1000;
end

%start to plot
if numel(t) < numel(Yr),
  warning('Protocol pulses defined require longer than the specified duration, cutting off at duration');
  Yr = Yr(1:numel(t));
  Yg = Yg(1:numel(t));
  Yb = Yb(1:numel(t));
end

% Package them up for return
result = [Yr Yg Yb] ;
