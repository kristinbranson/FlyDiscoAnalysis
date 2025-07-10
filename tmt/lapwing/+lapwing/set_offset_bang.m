function set_offset_bang(h, offset) 
    old_position = h.Position ;
    old_size = old_position(3:4) ;
    new_position = [offset old_size] ;
    h.Position = new_position ;
end
