function [firstframes, firstframes_off, endframes_off, local_nframes, local_indicatorframes] = ...
    DetermineStartAndEndOfResultsMovieSnippets(commonregistrationparams, registration_params, ctraxresultsmovie_params, local_nframes_input, ...
                                               local_indicatorframes_input, defaultparams, indicatorstruct, ledprotocolstruct, metadata)    
  
  local_nframes = local_nframes_input ;
  local_indicatorframes = local_indicatorframes_input ;
  if ~commonregistrationparams.OptogeneticExp,
    nframes = registration_params.end_frame-registration_params.start_frame + 1;
    firstframes_off = min(max(0,round(ctraxresultsmovie_params.firstframes*nframes)),nframes-1);
    firstframes_off(ctraxresultsmovie_params.firstframes < 0) = nan;
    middleframes_off = round(ctraxresultsmovie_params.middleframes*nframes);
    middleframes_off(ctraxresultsmovie_params.middleframes < 0) = nan;
    endframes_off = round(ctraxresultsmovie_params.endframes*nframes);
    endframes_off(ctraxresultsmovie_params.endframes < 0) = nan;
    idx = ~isnan(middleframes_off);
    firstframes_off(idx) = ...
     min(nframes-1,max(0,middleframes_off(idx) - ceil(local_nframes(idx)/2)));
    idx = ~isnan(endframes_off);
    firstframes_off(idx) = ...
     min(nframes-1,max(0,endframes_off(idx) - local_nframes(idx)));
    endframes_off = firstframes_off + local_nframes - 1;
    firstframes = registration_params.start_frame + firstframes_off;
  else
    if defaultparams,
      indicatorLED = indicatorstruct.indicatorLED ;
      protocol = ledprotocolstruct.protocol ;
      if strcmp(metadata.assay,'FlyBubbleRGB') || strcmp(metadata.assay,'FlyBowlRGB')
          if isfield(protocol,'Rintensity')
              RGBprotocol = protocol;
              clear protocol;
              % test if RGBprotocol has only one active color
              countactiveLEDs = [double(any(RGBprotocol.Rintensity));double(any(RGBprotocol.Gintensity));double(any(RGBprotocol.Bintensity))];
              % check that there is 1 and only 1 color LED used in protocol
              if sum(countactiveLEDs) == 0
                  error('ChR = 1 for LED protcol with no active LEDs')
              elseif sum(countactiveLEDs) > 1
                  error('More than one active LED color in protocol. Not currently supported')
              end
              % call function that transforms new protocol to old protocol
              protocol = ConvertRGBprotocol2protocolformat(RGBprotocol,countactiveLEDs) ;
          end
      end

      nsteps = numel(protocol.stepNum);
      % if there's only 1 step 
      if nsteps == 1,
        iter = protocol.iteration;
        n = ceil(iter/6);     % take every nth iteration of the stimulus
        indicatorframes = 1:n:iter;
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

      local_indicatorframes = indicatorframes;
      firstframes_off = indicatorLED.startframe(indicatorframes) - ctraxresultsmovie_params.nframes_beforeindicator;  
      % need to do this for specific files as well
      local_nframes = ones(1,length(firstframes_off))*local_nframes;
      endframes_off = firstframes_off + local_nframes -1;
      firstframes = registration_params.start_frame + firstframes_off;
    else
      if ~isempty(indicatorstruct) , 
        indicatorLED = indicatorstruct.indicatorLED ;
        is_local_indicatorframes_valid = (local_indicatorframes <= length(indicatorLED.startframe)) ;
        if ~all(is_local_indicatorframes_valid) ,
            error(['The largest value in local_indicatorframes (%d) is too big for indicatorLED.startframe (length = %d).  ' ... 
                   'Likely there was a problem with indicator detection.'], ...
                    max(local_indicatorframes), length(indicatorLED.startframe)) ;
            local_indicatorframes = local_indicatorframes(is_local_indicatorframes_valid) ;   
        end
        firstframes_off = indicatorLED.startframe(local_indicatorframes) - ctraxresultsmovie_params.nframes_beforeindicator;
        endframes_off = firstframes_off + local_nframes -1 ;
        firstframes = registration_params.start_frame + firstframes_off;
        local_nframes = ones(1,length(firstframes_off)).*local_nframes;
      end
    end
  end
end
