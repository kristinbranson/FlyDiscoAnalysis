function result = struct_from_name_value_list(name_value_list)
  % Convert a (flattened) list of name, value pairs to a struct, with the names
  % becoming field names.
  
  % Parse the default name-value list
  names = name_value_list(1:2:end) ;
  values = name_value_list(2:2:end) ;
  if length(names) ~= length(values) ,
    error('Odd number of elements in the name-value list') ;
  end  
  if ~all(cellfun(@ischar, names)) ,
    error('All names in the name-value list must be char arrays') ;
  end    
  
  field_count = length(names) ;
  result = struct() ;
  for i = 1 : field_count ,
    name = names{i} ;
    value = values{i} ;
    result.(name) = value ;
  end
end
