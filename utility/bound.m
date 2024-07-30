function result = bound(x, x_min, x_max)
    result = max(x_min, min(x, x_max, 'includenan'), 'includenan') ;
end
