function ctrax_results = ctrax_results_struct_from_flytracker_jaaba_trx_struct(trx_file_as_struct)    
    % Ctrax file should look like this:    
    %   startframe          1x1                  8  int64               
    %   timestamps       3624x1              28992  double              
    %   ntargets         3624x1              28992  double              
    %   identity        54361x1             434888  double              
    %   x_pos           54361x1             434888  double              
    %   y_pos           54361x1             434888  double              
    %   angle           54361x1             434888  double              
    %   maj_ax          54361x1             434888  double              
    %   min_ax          54361x1             434888  double              
    %
    % In the ctrax output, timestamps and ntargets are frame_count x 1
    % Other things have one element per (frame,individual) pair, for all
    % the indivuduals that could be detected in each frame.
    % ntargets gives the number of individuals detected in each frame.
    % Thus sum(ntargets)==pair_count.
    % identity, x_pos, y_pos, angle, maj_ax, and min_ax are are pair_count
    % x 1.  identity gives the individual identifier (an integer) for each pair.
    % x_pos, y_pos are in pixels.  angle is in radians.  maj_ax, min_ax are
    % in pixels.

    %trx_file_as_struct = load(flytracker_jaaba_trx_file_name, '-mat') ;
    trx = trx_file_as_struct.trx ;
    timestamps = trx_file_as_struct.timestamps(:) ;  % want col vector
    clear('trx_file_as_struct') ;  % no longer needed
    fly_count = length(trx) ;
    frame_count = length(timestamps) ;

    % Set the time of the first frame to 0, b/c that's how Ctrax rolls
    if frame_count>0 ,
        timestamps = timestamps-timestamps(1) ;
    end
    
    % Create the structure to populate
    ctrax_results = struct() ;
    
    % Not sure what startframe is supposed to be
    startframe = 0 ;
    ctrax_results.startframe = startframe ;  % is this right?

    % Store the timestamps
    ctrax_results.timestamps = timestamps ;
    
    % Make a big matrix indicating which flies are present in each frame,
    % and also the major variables we want 
    is_fly_i_present_in_frame_j = false(fly_count, frame_count) ;
    x_for_fly_i_from_frame_j = nan(fly_count, frame_count) ;
    y_for_fly_i_from_frame_j = nan(fly_count, frame_count) ;
    theta_for_fly_i_from_frame_j = nan(fly_count, frame_count) ;
    a_for_fly_i_from_frame_j = nan(fly_count, frame_count) ;
    b_for_fly_i_from_frame_j = nan(fly_count, frame_count) ;
    for fly_index = 1 : fly_count, 
        firstframe = trx(fly_index).firstframe ;
        endframe = trx(fly_index).endframe ;
        is_fly_i_present_in_frame_j(fly_index, firstframe:endframe) = true ;
        x_for_fly_i_from_frame_j(fly_index, firstframe:endframe) = trx(fly_index).x ;
        y_for_fly_i_from_frame_j(fly_index, firstframe:endframe) = trx(fly_index).y ;
        theta_for_fly_i_from_frame_j(fly_index, firstframe:endframe) = trx(fly_index).theta ;
        a_for_fly_i_from_frame_j(fly_index, firstframe:endframe) = trx(fly_index).a ;
        b_for_fly_i_from_frame_j(fly_index, firstframe:endframe) = trx(fly_index).b ;
    end
    
    % Determine how many flies present in each frame
    ntargets = sum(is_fly_i_present_in_frame_j, 1)' ;
    ctrax_results.ntargets = ntargets ;  

    % Figure out how pairs correspond to fly indices and frame indices    
    % All the pairs for frame 1 come before those for frame 2, etc
    pair_count = sum(sum(is_fly_i_present_in_frame_j)) ;
    fly_index_for_fly_i_from_frame_j = repmat( (1:fly_count)', [1 frame_count]) ;
    fly_index_from_pair_index = fly_index_for_fly_i_from_frame_j(is_fly_i_present_in_frame_j) ;   
    
    % Determine the fly id code for each pair
    fly_id_code_from_fly_index = [trx(:).id]' ;    
    fly_id_code_from_pair_index = fly_id_code_from_fly_index(fly_index_from_pair_index) ;
    ctrax_results.identity = fly_id_code_from_pair_index ;  

    % Get the x,y,theta,a,b data out of the matrices in "frame-major" order
    x_from_pair_index = x_for_fly_i_from_frame_j(is_fly_i_present_in_frame_j) ;
    y_from_pair_index = y_for_fly_i_from_frame_j(is_fly_i_present_in_frame_j) ;
    theta_from_pair_index = theta_for_fly_i_from_frame_j(is_fly_i_present_in_frame_j) ;
    a_from_pair_index = a_for_fly_i_from_frame_j(is_fly_i_present_in_frame_j) ;
    b_from_pair_index = b_for_fly_i_from_frame_j(is_fly_i_present_in_frame_j) ;
    assert(length(x_from_pair_index) == pair_count) ;

    % Stuff the ellipse parameters into ctrax_results
    ctrax_results.x_pos = x_from_pair_index ;
    ctrax_results.y_pos = y_from_pair_index ;
    ctrax_results.angle = theta_from_pair_index ;
    ctrax_results.maj_ax = a_from_pair_index ;
    ctrax_results.min_ax = b_from_pair_index ;
end
