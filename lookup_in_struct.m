function result = lookup_in_struct(default_parameters, probe_name, fallback_value)
  % We assume default_parameters is a struct of parameters.  We look up the
  % probe name in the struct.  If not found, we return the fallback value,
  % unless no fallback value was given, in which case we error.
  
  try 
    result = default_parameters.(probe_name) ;
  catch me 
    if isequal(me.identifier, 'MATLAB:nonExistentField') ,
      if exist('fallback_value', 'var') ,
          result = fallback_value ;
      else
          error('lookup_in_struct:not_found', 'Probe name "%s" not found in struct, and no fallback value given', probe_name) ;
      end      
    else
      rethrow(me) ;
    end
  end
  
end
