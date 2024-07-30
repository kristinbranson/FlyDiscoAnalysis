function s = merge_structs(varargin)
  % s = merge_structs(s1,s2,s3,...)
  % Merge structures into one combined structure.
  % Rightmost argument has priority if both share a field name.
  % This one doesn't warn about common field names.
  
  s = struct();
  for argument_index = 1 : length(varargin) ,
    this_struct = varargin{argument_index} ;
    assert(isscalar(this_struct) && isstruct(this_struct), ...
           'All input arguments must be scalar structures.');
    field_names = fieldnames(this_struct) ; 
    field_count = length(field_names) ;
    for field_index = 1 : field_count ,
      field_name = field_names{field_index} ;
      s.(field_name) = this_struct.(field_name) ;
    end
  end
end
