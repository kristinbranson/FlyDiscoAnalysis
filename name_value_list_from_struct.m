function result = name_value_list_from_struct(s)
  % Convert a struct to a (flattened) list of (name, value) pairs.
  % Result will be a cell array with twice at many elements as s has fields.
  
  % Parse the default name-value list
  names = fieldnames(s) ;
  name_count = length(names) ;
  item_count = 2 * name_count ;
  result = cell(1, item_count) ;
  for name_index = 1 : name_count ,
    name = names{name_index} ;
    value = s.(name) ;
    result{2*name_index-1} = name ;
    result{2*name_index  } = value ;
  end
end
