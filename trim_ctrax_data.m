function trim_ctrax_data(output_mat_name, input_mat_name, frame_count)
%   startframe          1x1                  8  int64               
%   timestamps       3624x1              28992  double              
%   ntargets         3624x1              28992  double              
%   identity        54361x1             434888  double              
%   x_pos           54361x1             434888  double              
%   y_pos           54361x1             434888  double              
%   angle           54361x1             434888  double              
%   maj_ax          54361x1             434888  double              
%   min_ax          54361x1             434888  double              

    in = load(input_mat_name) ;
    startframe_in = in.startframe ;  % scalar
    timestamps_in = in.timestamps ;  % frame_count_in x 1 
    ntargets_in = in.ntargets ;  % frame_count_in x 1   
    identity_in = in.identity ;  % pair_count_in x 1
    x_pos_in = in.x_pos ;  % pair_count_in x 1
    y_pos_in = in.y_pos ;  % pair_count_in x 1
    angle_in = in.angle ;  % pair_count_in x 1    
    maj_ax_in = in.maj_ax ;  % pair_count_in x 1
    min_ax_in = in.min_ax ;  % pair_count_in x 1
    
    framenumber = construct_frame_index_from_pair_index(ntargets_in) ;
    is_keeper_from_pair_index = (framenumber<=frame_count) ;
    
    startframe = startframe_in ;
    timestamps = timestamps_in(1:frame_count) ;
    ntargets = ntargets_in(1:frame_count) ;
    identity = identity_in(is_keeper_from_pair_index) ;
    x_pos = x_pos_in(is_keeper_from_pair_index) ;
    y_pos = y_pos_in(is_keeper_from_pair_index) ;
    angle = angle_in(is_keeper_from_pair_index) ;
    maj_ax = maj_ax_in(is_keeper_from_pair_index) ;
    min_ax = min_ax_in(is_keeper_from_pair_index) ;
    
    save(output_mat_name, 'startframe', 'timestamps', 'ntargets', 'identity', 'x_pos', 'y_pos', 'angle', 'maj_ax', 'min_ax') ;
end



function result = construct_frame_index_from_pair_index(target_count_from_frame_index)
    % For each frame, there is some number of individuals tracked.  For
    % frame i, this number is given by target_count_from_frame_index(i).
    % In the main ctrax output arrays, like x_pos, etc, all the <frame, id>
    % pairs for the first frame come first, then for the second frame, etc.
    % And for frame i, there are target_count_from_frame_index(i) pairs.
    % So we can build up an array that represents the frame index given the
    % pair index.  This is the frame_index_from_pair_index array, which is
    % our result here.
    
    pair_count = sum(target_count_from_frame_index) ;
    result = zeros(1, pair_count);
    pair_count_so_far = 0;
    frame_count = length(target_count_from_frame_index) ;
    for frame_index = 1:frame_count ,
        target_count = target_count_from_frame_index(frame_index) ;
        result(pair_count_so_far+1:pair_count_so_far+target_count) = frame_index;
        pair_count_so_far = pair_count_so_far + target_count ;
    end
    
    assert(pair_count_so_far==pair_count) ;  % sanity check
    % result(pair_index) now gives the (one-based) frame index for the <frame, id>
    % pair with index pair_index.    
end
