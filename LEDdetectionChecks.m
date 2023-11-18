function [passed,errormsg] = LEDdetectionChecks(protocol,indicatordata)

% check if number of LED bouts detected in video matches protocol
% for RGB protocols
% outputs
% passed = 1 logical for single color
% passed - 1x3 logical for RGB

if isfield(protocol,'Rintensity')
    % three color protocol
    countactiveLEDs = [double(any(protocol.Rintensity));double(any(protocol.Gintensity));double(any(protocol.Bintensity))];
    
    colorindex = 1:3;
    colors = {'Red','Green','Blue'};
    colors_iterations = {'Riteration','Giteration','Biteration'};
    passed = nan(1,3);
    for j = 1:numel(colorindex(logical(countactiveLEDs)))
        idx = colorindex(logical(countactiveLEDs));
        i = idx(j);
        totalprotocolbouts = sum(protocol.(colors_iterations{i}));
        if numel(indicatordata.startframe) ~= numel(indicatordata.endframe)
            
        end

        if totalprotocolbouts == numel(indicatordata.startframe)
            passed(i) = true;
            errormsg{i} = {};
        else
            passed(i) = false;
            errormsg{i} = sprintf('Num of %s LED stims: %d do not match protocol: %d',colors{i}, numel(indicatordata.startframe), totalprotocolbouts);
        end

    end
% for single color protocol
else
    % single color protocol
    totalprotocolbouts = sum(protocol.iteration);
    assert(numel(indicatordata.startframe) == numel(indicatordata.endframe),'movie starts or ends with indicator LED on')
    if totalprotocolbouts ~= numel(indicatordata.startframe)
        passed = false;
        errormsg = sprintf('Num of indicator stims: %d does not match protocol: %d', numel(indicatordata.startframe), totalprotocolbouts);
    else
        passed = true;
        errormsg = {};
    end
end