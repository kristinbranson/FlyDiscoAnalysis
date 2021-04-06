function result = short_host_name_from_host_name(host_name)
    parts = strsplit(host_name, '.') ;
    if isempty(parts) ,
        error('Unable to derive a short host name from host name "%s"'. host_name) ;
    else
        first_part = parts{1} ;
        if isempty(first_part) ,
            error('Unable to derive a short host name from host name "%s"', host_name) ;
        else
            result = first_part ;
        end
    end            
end
