function result = is_fid_valid(fid)
    result = ~isempty(fopen(fid)) ;
end
