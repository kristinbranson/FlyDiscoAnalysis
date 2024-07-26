function result = ensafen_readframe_fcn(readframe1, nframes, nr, nc, dataclass)
    % Take an unsafe frame-reading function and return a safer one.
    % The returned function returns black frames if an out-of-range 
    % frame index is requested.
    function frame = readframe2(frame_index)
        if 1<=frame_index && frame_index<nframes ,
            frame = feval(readframe1, frame_index) ;
        else
            frame = zeros([nr nc], dataclass) ;
        end
    end
    result = @readframe2 ;
end
