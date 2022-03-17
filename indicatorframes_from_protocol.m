function indicatorframes = indicatorframes_from_protocol(protocol)
    % Given the protocol, determine which stimuli to use for each snippet.
    % indicatorframes(snippet_index) gives the stimulus index for that snippet
    nsteps = numel(protocol.stepNum);
    % if there's only 1 step
    if nsteps == 1,
        iter = protocol.iteration;
        n = ceil(iter/6);     % take every nth iteration of the stimulus
        indicatorframes = 1:n:iter;  % "indicatorframes" could maybe also be named step_index_from_snippet_index, I think  -- ALT, 2022-03-16
    elseif nsteps <= 3,
        % assumes steps have more 2 or more iterations
        if ~all([protocol.iteration] > 1)
            error('Step with iteration < 2. User needs to make a specific ctrax results movie param file')
        end
        indicatorframes = zeros(1,2*nsteps);
        for step = 1:nsteps,
            iter = protocol.iteration(step);
            if step==1,
                indicatorframes(1)=1;
                indicatorframes(2)=iter;
                %         this logic doesn't work properly
                %         elseif
                %           indicatorframes(2*step-1)=indicatorframes(2*step-2)+1;
                %           indicatorframes(2*step)=indicatorframes(2*step-1);
            elseif step == 2
                indicatorframes(3) = indicatorframes(2)+1;
                indicatorframes(4) = indicatorframes(2)+iter;
            elseif step == 3
                indicatorframes(5) = indicatorframes(4)+1;
                indicatorframes(6) = indicatorframes(4)+iter;

            end
        end
    else
        n = ceil(nsteps/6);
        index = 1;
        indicatorframes=ones(1,numel(1:n:nsteps));
        for step = 1:n:nsteps,
            if step==1,
                indicatorframes(index)=1;
                index=index+1;
            else
                indicatorframes(index) = indicatorframes(index-1);
                for i = step-n:step-1
                    indicatorframes(index) = indicatorframes(index)+protocol.iteration(i);
                end
                index=index+1;

            end
        end
    end

    % make sure none of the steps are blank (have no stimulus)
    stim = 0;
    step = 0;
    for i=1:numel(indicatorframes),
        while indicatorframes(i) > stim,
            step=step+1;
            stim = stim+protocol.iteration(step);
        end
        if (protocol.intensity(step) == 0),
            if ~(i==numel(indicatorframes)),
                error('Step with intensity = 0 in middle of experiment. User needs to make specific ctrax results movie params file')
            end
            indicatorframes(i) = [];
        end
    end
    
end
