function result = lookup_in_name_value_list(default_parameters, probe_name, fallback_value)
  % We assume default_parameters in a (row) list of name, value
  % pairs.  We look up the probe name in the list.  If not found, we throw an
  % error.
  
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
  
  % Find the probe_name
  i = find(strcmp(probe_name, default_names), 1) ;
  if isempty(i) ,
      if exist('fallback_value', 'var') ,
          result = fallback_value ;
      else
          error('lookup_in_name_value_list:not_found', 'Probe name "%s" not found in name-value list', probe_name) ;
      end
  else      
      % Return the corresponding value
      result = default_values{i} ;
  end
end
