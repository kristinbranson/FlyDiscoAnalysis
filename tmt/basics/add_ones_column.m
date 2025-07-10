function result = add_ones_column(xyz)
    m = size(xyz, 1) ;
    result = [xyz ones(m,1)] ;
end
