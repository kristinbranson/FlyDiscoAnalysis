function [firstframes, firstframes_off, endframes_off, nframes, indicatorframes] = ...
        DetermineStartAndEndOfResultsMovieSnippets(isOptogeneticExp, registration_data, ctraxresultsmovie_params, ...
                                                   is_using_default_ctrax_results_movie_params, indicator_data, raw_protocol, ...
                                                   indicator_params)

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
    %          snippet for optogenetic experiments.  (Not the frame index, the index
    %          of the stimulus itself within the sequence of stimuli.)  For
    %          non-optogenetic experiments, always set to [].
    
    nframes_from_params = ctraxresultsmovie_params.nframes ;
    
    if isOptogeneticExp ,
        % If an optogenetic experiment
        if is_using_default_ctrax_results_movie_params,
            % Extract the indicatorframes from the protocol
            protocol = downmixProtocolIfNeeded(raw_protocol, indicator_params) ;            
            indicatorframes = indicatorframes_from_protocol(protocol) ;            
        else
            % Get the indicatorframes from the ctrax-results params
            indicatorframes = ctraxresultsmovie_params.indicatorframes ;
        end
        % "indicatorframes" could maybe also be named train_index_from_snippet_index,
        % if one was being long-winded, I think  -- ALT, 2024-03-22
        indicatorLED = downmix_indicatorLED(indicator_data.indicatorLED) ;
        startframe = indicatorLED.startframe ;
        indicator_train_count =  numel(startframe) ;
        are_all_indicatorframes_valid = (max(indicatorframes) <= indicator_train_count) ;
        if ~all(are_all_indicatorframes_valid) ,
            error(['The largest train index in the ctrax-results movie params (%d) is larger than the number of trains in the LED indicator trace (%d).  ' ...
                   'Likely there was a problem with indicator detection.'], ...
                  max(indicatorframes), ...
                  indicator_train_count) ;
        end
        firstframes_off = startframe(indicatorframes) - ctraxresultsmovie_params.nframes_beforeindicator;        
        % Promote nframes to a vector if necessary
        if isvector(firstframes_off) && isscalar(nframes_from_params) ,
            nframes = repmat(nframes_from_params, size(firstframes_off)) ;
        else
            nframes = nframes_from_params ;
        end
        endframes_off = firstframes_off + nframes - 1 ;
    else
        % If not an optogenetic experiment        
        % In this case, nframes_from_params should be the same size as
        % ctraxresultsmovie_params.firstframes, and their common length is the snippet
        % count.        
        tracked_frame_count = registration_data.end_frame - registration_data.start_frame + 1;
        % We define the "normalized position" within the tracked frames to be a number
        % on the interval [0,1], where 0 means the first frame, and 1 means the last
        % frame, and 0.5 means the middle frame.  You get the idea.
        firstframes_from_params = ctraxresultsmovie_params.firstframes ;  % e.g. [ 0 -1 -1 ]
          % firstframes_from_params(i) is the normalized position of snippet i, if snippet
          % i is to be a start-anchored snippet.  Otherwise, firstframes_from_params(i) is
          % -1.  A "start-anchored" snippet is one where the given normalized position
          % specifies where the snippet should *start*.
        middleframes_from_params = ctraxresultsmovie_params.middleframes ;  % e.g. [ -1 0.5 -1 ] 
          % middleframes_from_params(i) is the normalized position of snippet i, if snippet
          % i is to be a middle-anchored snippet.  Otherwise, middleframes_from_params(i) is
          % -1.  A "middle-anchored" snippet is one where the given normalized position
          % specifies where the middle of the snippet should be.
        endframes_from_params = ctraxresultsmovie_params.endframes ;  % e.g. [ -1 -1 1 ] 
          % endframes_from_params(i) is the normalized position of snippet i, if snippet
          % i is to be a end-anchored snippet.  Otherwise, endframes_from_params(i) is
          % -1.  An "end-anchored" snippet is one where the given normalized position
          % specifies where the snippet should *end*.        
        unbounded_firstframes_off = ...
            fifelse(firstframes_from_params >= 0, round(endframes_from_params*tracked_frame_count), ...
                    middleframes_from_params >= 0, round(middleframes_from_params*tracked_frame_count - nframes_from_params/2), ...
                    endframes_from_params >= 0, round(endframes_from_params*tracked_frame_count - nframes_from_params), ...
                    nan(size(firstframes_from_params))) ;
        firstframes_off = bound(unbounded_firstframes_off, 0, tracked_frame_count-1) ;             
        endframes_off = firstframes_off + nframes_from_params - 1 ;
        indicatorframes = [] ;  
            % indicatorframes is not used for non-optogenetic experiments, so doesn't really
            % matter what we return here.
        nframes = nframes_from_params ;
    end
    firstframes = registration_data.start_frame + firstframes_off;
end  % function
