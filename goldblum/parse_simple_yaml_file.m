function result = parse_simple_yaml_file(file_name)
    % Parse a file that is a list of field_name : value pairs, one per line.
    % Return the result in a flattened cell array of field_name, value pairs.  
    
    lines = read_file_into_cellstring(file_name) ;
    result = cell(1,0) ;
    for i = 1 : length(lines) ,
        raw_line = lines{i} ;
        line  = strtrim(raw_line) ;
        if isempty(line) ,
            continue
        end        
        colon_indices = strfind(line, ':') ;
        if isempty(colon_indices) , 
            error('At line %d: Unable to parse "%s"', i, raw_line) ;
        end
        colon_index = colon_indices(1) ;
        raw_field_name = line(1:(colon_index-1)) ;
        if length(line) < colon_index+1 ,
            error('At line %d: Unable to parse "%s"', i, raw_line) ;
        end            
        raw_value_as_string = line((colon_index+1):end) ;
        field_name = strtrim(raw_field_name) ;
        field_value_as_string = strtrim(raw_value_as_string) ;
        value = eval(field_value_as_string) ;  % risky
        result = horzcat(result, {field_name, value}) ;  %#ok<AGROW>
    end
end
