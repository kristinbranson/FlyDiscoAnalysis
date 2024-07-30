function result = str2logical(str) 
    trimmed_str = strtrim(str) ;
    if strcmp(trimmed_str, 'true') ,
        result = true ;
    elseif strcmp(trimmed_str, 'false') ,
        result = false ;
    else
        result = logical(str2double(str)) ;
    end        
end
