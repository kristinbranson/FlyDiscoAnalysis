function result = coerce_to_on_off_force(whatever)
    % Coerce lots of different things to one of 'on', 'off', 'force'
    if isstring(whatever) ,
        whatever = char(whatever) ;
    end
    if ischar(whatever) ,
        lowered_trimmed_string = lower(strtrim(whatever)) ;
        if strcmp(lowered_trimmed_string, 'on') || strcmp(lowered_trimmed_string, 'true') ,
            result = 'on' ;
        elseif strcmp(lowered_trimmed_string, 'off')  || strcmp(lowered_trimmed_string, 'false') ,
            result = 'off' ;
        elseif strcmp(lowered_trimmed_string, 'force') ,
            result = 'force' ;
        else
            error('Don''t know how to convert char array ''%s'' to on/off/force', whatever) ;
        end        
    elseif islogical(whatever) ,
        if whatever ,
            result = 'on' ;
        else
            result = 'off' ;
        end
    elseif isnumeric(whatever) ,
        whatever = logical(whatever) ;
        if whatever ,
            result = 'on' ;
        else
            result = 'off' ;
        end                
    else        
        error('Don''t know how to coerce something to on/off/force') ;
    end        
end
