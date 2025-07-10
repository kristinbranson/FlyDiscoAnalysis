function result = add_fields(s, varargin) 
    % A pure-functional way to add fields to a scalar struct.
    result = s ;
    vararg_count = length(varargin) ;
    pair_count = round(vararg_count/2) ;
    if vararg_count ~= 2*pair_count ,
        error('Number of arguments to %s must be odd', mfilename()) ;
    end
    for pair_index = 1 : pair_count ,
        field_name_vararg_index = 2*(pair_index-1) + 1 ;
        field_value_vararg_index = 2*(pair_index-1) + 2 ;
        field_name = varargin{field_name_vararg_index} ;
        field_value = varargin{field_value_vararg_index} ;
        result.(field_name) = field_value ;
    end
end
