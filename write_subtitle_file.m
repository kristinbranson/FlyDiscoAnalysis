function subtitlefile = ...
        write_subtitle_file(expdir, nframes, isOptogeneticExp, ctraxresultsmovie_params, basename, firstframes_off, endframes_off, ...
                            raw_protocol, indicatorstruct, indicatorframes)
    
    snippet_count = length(nframes) ;
    
    subtitlefile = fullfile(expdir,'subtitles.srt');
    if exist(subtitlefile,'file'),
        delete(subtitlefile);
    end
    fid = fopen(subtitlefile,'w');
    if fid < 0 ,
        error('Unable to open file %s for writing', subtitlefile) ;
    end
    cleaner = onCleanup(@()(fclose(fid))) ;
    
    dt = [0,nframes];
    ts = cumsum(dt);
    
    if isOptogeneticExp ,
        protocol = downmixProtocolIfNeeded(raw_protocol) ;
        
        indicatorLED = indicatorstruct.indicatorLED ;
        stimtimes = indicatorLED.starttimes(indicatorframes);
        stimulus_count = length(indicatorLED.starttimes) ;

        % Determine the step index for each snippet
        step = zeros(1,snippet_count) ;
        step_index = 1 ;  
        t = protocol.duration(step_index)/1000;        
        for snippet_index = 1:snippet_count,
            while stimtimes(snippet_index) > t
                step_index = step_index+1;
                t = t + protocol.duration(step_index)/1000;
            end            
            step(snippet_index) = step_index;
        end
            
        % For each snippet, compute several stimulus properties
        freq = zeros(1,snippet_count) ;
        intensity = zeros(1,snippet_count) ;
        dutycycle = zeros(1,snippet_count) ;
        duration = zeros(1,snippet_count) ;
        stim_type = cell(1, snippet_count) ;
        for snippet_index = 1:snippet_count,
            step_index = step(snippet_index) ;
            
            % not sure this logic is good
            if protocol.pulseNum(step_index) > 1,
                this_freq = 1/(protocol.pulsePeriodSP(step_index)/1000) ;
                freq(snippet_index) = this_freq ;
                stim_type{snippet_index} = ['Plsd ', num2str(this_freq), 'Hz'];
            else
                freq(snippet_index) = 0 ;                
                stim_type{snippet_index} = 'Cnst';
            end
            
            intensity(snippet_index) = protocol.intensity(step_index) ;
            dutycycle(snippet_index) = (protocol.pulseWidthSP(step_index)/protocol.pulsePeriodSP(step_index))*100 ;
            duration(snippet_index) = protocol.pulseNum(step_index)*protocol.pulsePeriodSP(step_index) ;            
        end
        
        for snippet_index = 1:snippet_count ,
            fprintf(fid,'%d\n',snippet_index);
            
            t_start = ts(snippet_index)/ctraxresultsmovie_params.fps/(3600*24);
            t_end = (ts(snippet_index) + (ts(snippet_index+1)-1-ts(snippet_index))/ctraxresultsmovie_params.subdecimationfactor)/ ...
                ctraxresultsmovie_params.fps/(3600*24);
            
            fprintf(fid,'%s --> %s\n',...
                datestr(t_start,'HH:MM:SS,FFF'),...
                datestr(t_end,'HH:MM:SS,FFF'));
            fprintf(fid,'%s\n%s %s %s %s %s %s\n\n',basename,...
                ['Step ', num2str(step(snippet_index))],...
                ['(Stim ', num2str(indicatorframes(snippet_index)),'/', ...
                num2str(stimulus_count),'):'],...
                stim_type{snippet_index},...
                ['(',num2str(dutycycle(snippet_index)),'% on)'],...
                ['at ',num2str(intensity(snippet_index)), '% intensity'],...
                ['for ',num2str(duration(snippet_index)), ' ms']);
        end
    else
        % non-optogenetic experiment
        for snippet_index = 1:snippet_count ,
            fprintf(fid,'%d\n',snippet_index);
            fprintf(fid,'%s --> %s\n',...
                datestr(ts(snippet_index)/ctraxresultsmovie_params.fps/(3600*24),'HH:MM:SS,FFF'),...
                datestr((ts(snippet_index+1)-1)/ctraxresultsmovie_params.fps/(3600*24),'HH:MM:SS,FFF'));
            fprintf(fid,'%s, fr %d-%d\n\n',basename,...
                firstframes_off(snippet_index)+1,...
                endframes_off(snippet_index)+1);
        end
    end
end
    
