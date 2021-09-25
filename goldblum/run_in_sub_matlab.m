function stdout = run_in_sub_matlab(do_actually_shell_out, function_handle, varargin)
    % Call a matlab function by invoking a sub-matlab at the command line to execute
    % it.  Doesn't provide a way to get function return values, so only useful for
    % functions with side-effects.
    if do_actually_shell_out ,
        function_name = func2str(function_handle) ;
        arg_string = generate_arg_string(varargin{:}) ;
        matlab_command = sprintf('modpath; %s(%s);', function_name, arg_string) ;
        %matlab_command = sprintf('%s(%s);', function_name, arg_string) ;
        bash_command_as_list = {'matlab', '-batch', matlab_command} ;
        stdout = system_from_list_with_error_handling(bash_command_as_list) ;
    else
        % Just call the function normally
        feval(function_handle, varargin{:}) ;
    end        
end



function result = tostring(thing)
    % Converts a range of things to strings that will eval to the thing
    if ischar(thing) ,
        result = sprintf('''%s''', thing) ;
    elseif isnumeric(thing) || islogical(thing) ,
        result = mat2str(thing) ;
    elseif isstruct(thing) && isscalar(thing) ,
        result = 'struct(' ;
        field_names = fieldnames(thing) ;
        field_count = length(field_names) ;
        for i = 1 : field_count ,
            field_name = field_names{i} ;
            field_value = thing.(field_name) ;
            field_value_as_string = tostring(field_value) ;
            subresult = sprintf('''%s'', {%s}', field_name, field_value_as_string) ;
            result = horzcat(result, subresult) ; %#ok<AGROW>
            if i<field_count ,
                result = horzcat(result, ', ') ; %#ok<AGROW>
            end            
        end
        result = horzcat(result, ')') ;
    elseif iscell(thing) && (isempty(thing) || isvector(thing)) ,
        if isempty(thing) ,
            result = sprintf('cell(%d,%d)', size(thing,1), size(thing,2)) ;
        else
            if iscolumn(thing) ,
                separator = ';' ;
            else
                separator = ',' ;
            end            
            result = '{ ' ;
            element_count = length(thing) ;
            for i = 1 : element_count ,
                element_value = thing{i} ;
                element_value_as_string = tostring(element_value) ;
                result = horzcat(result, element_value_as_string) ; %#ok<AGROW>
                if i<element_count ,
                    result = horzcat(result, [separator ' ']) ; %#ok<AGROW>
                end            
            end
            result = horzcat(result, ' }') ;
        end
    else
        error('Don''t know how to convert something of class %s to string', class(thing)) ;
    end
end



function result = generate_arg_string(varargin) 
    arg_count = length(varargin) ;
    result = char(1,0) ;  % fall-through in case of zero args
    for i = 1 : arg_count ,
        this_arg = varargin{i} ;
        this_arg_as_string = tostring(this_arg) ;
        if i == 1 ,
            result = this_arg_as_string ;
        else
            result = horzcat(result, ', ', this_arg_as_string) ;  %#ok<AGROW>
        end
    end
end
