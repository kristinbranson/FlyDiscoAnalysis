function position_in_upper_left_corner_of_primary_screen_bang(fig, left_padding, top_padding)
    % What it says on the tin.  left_padding and top_padding should be in
    % pels.
    
    % get the screen size so we can position the figure window
    root_units=get(0,'Units');
    set(0,'Units','pixels');
    screen_dims=get(0,'ScreenSize');
    %screen_width=screen_dims(3);
    screen_height=screen_dims(4);
    set(0,'Units',root_units);
    
    % Make adjustments for platform
    if ismac() ,
        top_padding = top_padding + 22 ;  % Mac menu bar is 22 pels high
    end
    
    % get the figure size, in the OuterPosition sense
    position = get(fig,'OuterPosition') ;
    size = position(3:4) ;
    
    % Compute the new figure offset
    new_x = left_padding ;
    new_y = screen_height - size(2) - top_padding ;    
    new_offset = [new_x new_y] ;
    
    % Set the position, using the new offset with the original size
    new_position = [new_offset size] ;
    set(fig, 'OuterPosition', new_position) ;    
end