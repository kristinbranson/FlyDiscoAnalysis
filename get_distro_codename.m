function result = get_distro_codename()
    % Returns 'Ubuntu' if on Ubuntu, 'OracleLinux' if on Oracle Linux, 'Windows' if on Windows, and 'macOS' if on
    % macOS.
    if isunix() ,
        if exist('/usr/bin/lsb_release', 'file') ,
            [return_code, stdout] = system('/usr/bin/lsb_release -i -s') ;
            if return_code==0 ,
                raw_result = strtrim(stdout) ;
            else
                error('Unable to determine operating system---call to lsb_release failed') ;
            end
        else
            raw_result = get_name_from_etc_os_release() ;
        end
        result = massage_disto_name(raw_result) ;
    else
        if ispc() ,
            result = 'Windows' ;
        elseif ismac() ,
            result = 'macOS' ;
        else
            error('Unable to determine operating system') ;
        end
    end
end



function result = get_name_from_etc_os_release()
    % Gets the distro name from reading, parsing /etc/os-release.
    % Parses it in a half-assed way that nonetheless seems to be sufficient for
    % our needs.
    % Errors if /etc/os-release doesn't exist or we can't parse it, or the NAME
    % key is missing.
    if exist('/etc/os-release', 'file') ,
        line_from_line_number = read_file_into_cellstring('/etc/os-release') ;
        for line_number = 1 : length(line_from_line_number) ,
            line = line_from_line_number{line_number} ;
            maybe_kv_pair = parse_line(line) ;
            if isempty(maybe_kv_pair) 
                continue
            else
                kv_pair = maybe_kv_pair{1} ;
                key = kv_pair{1} ;
                if strcmp(key, 'NAME') ,
                    value = kv_pair{2} ;
                    % Strip quotes
                    name = strip(value, '"') ;
                    result = name ;
                    return
                else
                    continue
                end

            end
        end
        % If get here, failed to find NAME line 
        error('Unable to determine operating system---no lsb_release, no NAME line in /etc/os-release') ;                
    else
        error('Unable to determine operating system---no lsb_release, no /etc/os-release') ;        
    end
end



function result = parse_line(line)
    % Parse a single line of /etc/os-release, splitting at the first equals sign,
    % and trimming excess whitespace from the resulting strings.  Returns an
    % empty cell array if there is no equals sign.  Normally returns a singleton 
    % cell array whose only element is itself a two-element cell array containing key
    % and value.
    index_from_equals_index = strfind(line, '=') ;
    if isempty(index_from_equals_index) ,
        result = {} ;
    else
        index_of_first_equals = index_from_equals_index(1) ;
        raw_key = line(1:index_of_first_equals-1) ;
        if length(line) > index_of_first_equals ,
            raw_value = line(index_of_first_equals+1:end) ;
        else
            raw_value = '' ;
        end
        key = strtrim(raw_key) ;
        value = strtrim(raw_value) ;
        result = { { key, value } } ;
    end
end



function result = massage_disto_name(raw_name)
    % Different ways of getting distro name on Oracle Linux yield different
    % strings, so we sanitize here.
    if contains(raw_name, 'Oracle') ,
        result = 'OracleLinux' ;
    else
        result = raw_name ;
    end
end
