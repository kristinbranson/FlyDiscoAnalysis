function result = str2logical(raw_str)
    % Converts a string to a scalar logical value, without using eval()
    str = lower(strtrim(raw_str)) ;
    if strcmp(str, 'true') ,
        result = true ;
    elseif strcmp(str, 'false') ,
        result = false ;
    else
        % try converting to double, and accept any nonzero value as true
        as_double = str2double(str) ;
        if isreal(as_double) && ~isnan(as_double) ,  % isreal() test is so e.g. input 'i' errors rather than returning true
            result = (as_double~=0) ;
        else            
            % give up
            error('str2logical:bad_input', 'Can''t convert string "%s" to logical value', raw_str) ;
        end
    end
end
