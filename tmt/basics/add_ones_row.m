function result = add_ones_row(xyz)
    n = size(xyz, 2) ;
    result = [ xyz ; ...
               ones(1,n) ] ;
end
