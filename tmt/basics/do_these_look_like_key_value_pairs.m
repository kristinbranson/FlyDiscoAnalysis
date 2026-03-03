function result = do_these_look_like_key_value_pairs(args)
    % Test to see if args looks like a list of keyword-value pairs.
    % Returns true if args has an even number of elements, and all the
    % odd-index elements are strings.

    n = length(args) ;
    if mod(n, 2) ~= 0 ,
        result = false ;
        return
    end
    key_list = args(1:2:end) ;
    is_old_string_from_key_index = cellfun(@is_old_string, key_list) ;
    is_new_string_from_key_index = cellfun(@isstring, key_list) ;
    is_string_from_key_index = is_old_string_from_key_index | is_new_string_from_key_index ;
    result = all(is_string_from_key_index) ;
end


function result = is_old_string(x)
    result = ischar(x) && ismatrix(x) && size(x,1)<=1 ;
end
