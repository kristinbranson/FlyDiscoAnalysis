function [firstframes, firstframes_off, endframes_off, nframes, indicatorframes] = ...
        DetermineStartAndEndOfResultsMovieSnippets(commonregistrationparams, registration_params, ctraxresultsmovie_params, ...
                                                   input_nframes, input_indicatorframes, ...
                                                   is_using_default_ctrax_results_movie_params, indicator_file_contents, ledprotocol_file_contents, metadata)

    %DetermineStartAndEndOfResultsMovieSnippets  Use information from various
    %sources to determine which frames to use for each snippet included in the
    %results movie.  This function is a pure function.  In particular, it doesn't
    %read anything from disk.
    %
    %  On return:
    %      firstframes(snippet_index) gives the frame index of the first frame of the
    %          indicated snippet.
    %      firstframes_off(snippet_index) gives the frame *offset* of the first frame
    %          of the indicated snippet, relative to the start of the stimulus sequence in the
    %          video.  (The stimulus does not have to start at frame 1, and usually
    %          doesn't.)  The frame index of the start of the stimulus sequence must
    %          be provided in registration_params.start_frame.  Thus if the stimulus
    %          to be shown in snippet i had firstframes_off(i)==0, then
    %          firstframes(i) would be equal to registration_params.start_frame.
    %      endframes_off(snippet_index) gives the frame *offset* of the last frame of
    %          the indicated snippet.  Like with firstframes_off, these offsets are
    %          relative to registration_params.start_frame.
    %      nframes(snippet_index) gives the number of frames in the indicated
    %          snippet.  It will always be the case that 
    %          endframes_off == firstframes_off + nframes - 1.
    %      indicatorframes(snippet_index) gives the simulus index of the indicated
    %          snippet.  (Not the frame index, the index of the stimulus itself
    %          within the sequence of stimuli.)
    
    if commonregistrationparams.OptogeneticExp ,
        % if an optogenetic experiment
        if is_using_default_ctrax_results_movie_params,
            % Figure out the protocol given the metadata, LED protocol file contents
            protocol = determine_protocol(metadata, ledprotocol_file_contents) ;
            
            % Extract the indicatorframes from the protocol
            indicatorframes = indicatorframes_from_protocol(protocol) ;
            
            indicatorLED = indicator_file_contents.indicatorLED ;            
            firstframes_off = indicatorLED.startframe(indicatorframes) - ctraxresultsmovie_params.nframes_beforeindicator ;
        else
            % if do_use_default_parameters is false
            if isempty(indicator_file_contents) ,
                error('In %s(), do_use_default_parameters is false, but indicator_file_contents is empty.  This combination is not currently supported.', ...
                      mfilename()) ;
            else
                indicatorframes = input_indicatorframes ;
                indicatorLED = indicator_file_contents.indicatorLED ;
                is_local_indicatorframes_valid = (indicatorframes <= length(indicatorLED.startframe)) ;
                if ~all(is_local_indicatorframes_valid) ,
                    error(['The largest value in indicatorframes (%d) is too big for indicatorLED.startframe (length = %d).  ' ...
                        'Likely there was a problem with indicator detection.'], ...
                        max(indicatorframes), length(indicatorLED.startframe)) ;
                    indicatorframes = indicatorframes(is_local_indicatorframes_valid) ;
                end
                firstframes_off = indicatorLED.startframe(indicatorframes) - ctraxresultsmovie_params.nframes_beforeindicator;
            end
        end
        % Promote nframes to a vector if necessary
        if isvector(firstframes_off) && isscalar(input_nframes) ,
            nframes = repmat(input_nframes, size(firstframes_off)) ;
        else
            nframes = input_nframes ;
        end
        endframes_off = firstframes_off + nframes - 1 ;
    else
        % If not an optogenetic experiment
        tracked_frame_count = registration_params.end_frame-registration_params.start_frame + 1;
        firstframes_off = min(max(0,round(ctraxresultsmovie_params.firstframes*tracked_frame_count)),tracked_frame_count-1);
        firstframes_off(ctraxresultsmovie_params.firstframes < 0) = nan;
        middleframes_off = round(ctraxresultsmovie_params.middleframes*tracked_frame_count);
        middleframes_off(ctraxresultsmovie_params.middleframes < 0) = nan;
        endframes_off = round(ctraxresultsmovie_params.endframes*tracked_frame_count);
        endframes_off(ctraxresultsmovie_params.endframes < 0) = nan;
        idx = ~isnan(middleframes_off);
        firstframes_off(idx) = ...
            min(tracked_frame_count-1,max(0,middleframes_off(idx) - ceil(input_nframes(idx)/2)));
        idx = ~isnan(endframes_off);
        firstframes_off(idx) = ...
            min(tracked_frame_count-1,max(0,endframes_off(idx) - input_nframes(idx)));
        endframes_off = firstframes_off + input_nframes - 1;
        indicatorframes = input_indicatorframes ;  % TODO: Is this right?
        nframes = input_nframes ;
    end
    firstframes = registration_params.start_frame + firstframes_off;
end  % function
