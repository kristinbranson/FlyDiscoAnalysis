function result = char_array_from_value(value)
  % Returns a string representation of the value.
  % *Very* primitive right now.
  
  if ischar(value) && isrow(value) ,
    result = [ '''' value '''' ] ;
  elseif isstring(value) ,
    result = sprintf('"%s"', value) ;
  else
    if ~isscalar(value) ,
      error('%s only supports strings and scalar values', mfilename()) ;
    end
    if islogical(value) ,
      if value , 
        result = 'true' ;
      else
        result = 'false' ;
      end
    elseif isnumeric(value) , 
      result = num2str(value) ;
    else
      error('%s doesn''t support scalars of type', mfilename(), class(value)) ;
    end
  end
end
