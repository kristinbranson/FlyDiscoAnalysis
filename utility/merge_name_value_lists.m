function result = merge_name_value_lists(default_parameters, new_parameters)
  % We assume default_parameters and parameters are both (row) lists of name, value
  % pairs.  We merge them, with parameters taking precedence over
  % default_parameters.  The result is another list of name, value pairs.
  
  % Parse the default name-value list
  default_names = default_parameters(1:2:end) ;
  default_values = default_parameters(2:2:end) ;
  if length(default_names) ~= length(default_values) ,
    error('Odd number of elements in the default_parameters name-value list') ;
  end  
  if ~all(cellfun(@ischar, default_names)) ,
    error('All names in the default_parameters name-value list must be char arrays') ;
  end    
  %default_pair_count = length(default_values) ;
  
  % Parse the new name-value list
  new_names = new_parameters(1:2:end) ;
  new_values = new_parameters(2:2:end) ;  
  if length(new_names) ~= length(new_values) ,
    error('Odd number of elements in the new_parameters name-value list') ;
  end
  if ~all(cellfun(@ischar, new_names)) ,
    error('All names in the new_parameters name-value list must be char arrays') ;
  end  
  %new_pair_count = length(new_values) ;
  
  [result_names, ~, indices_of_default_names_in_result_as_col] = union(new_names, default_names, 'stable') ;
  result_pair_count = length(result_names) ;
  
  % Assemble the list of result values.  B/c we use stable sort, they will be all
  % the new_values, followed by a subset of the default_values.
  result_values = horzcat(new_values, default_values(indices_of_default_names_in_result_as_col') ) ;
  
  % Intercalate result_names and result_values
  result = cell(1, 2*result_pair_count) ;
  result(1:2:end) = result_names ;
  result(2:2:end) = result_values ;
end
