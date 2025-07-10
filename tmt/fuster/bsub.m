function job_id = bsub(do_actually_submit, ...
                       do_try, ...
                       slot_count, ...
                       stdouterr_file_name, ...
                       options, ...
                       do_use_xvfb, ...
                       submit_host_name, ...
                       function_handle, ...
                       varargin)
    % Wrapper for LSF bsub command.  Returns job id as a double.
    % Throws error if anything goes wrong.
    if isempty(do_actually_submit) ,
        do_actually_submit = true ;        
    end
    if isempty(do_try) ,
        do_try = true ;        
    end
    if isempty(slot_count) ,
        slot_count = 1 ;
    end    
    if isempty(stdouterr_file_name) ,
        stdouterr_file_name = '/dev/null' ;
    end    
    if isempty(do_use_xvfb) ,
        do_use_xvfb = false ;
    end
    if isempty(submit_host_name) ,
        submit_host_name = '' ;
    end
    if do_actually_submit ,
        function_name = func2str(function_handle) ;
        arg_string = generate_arg_string(varargin{:}) ;        
        matlab_command = sprintf('cd(''%s''); modpath; %s(%s);', pwd(), function_name, arg_string) ;
        matlab_version_as_string = matlab_version_string() ;  % e.g. 2022a
        matlab_executable_path = sprintf('/misc/local/matlab-%s/bin/matlab', matlab_version_as_string) ;
          % Use same version of Matlab as is being used to run this code.
        if do_use_xvfb ,
            bash_command = sprintf('/usr/bin/xvfb-run -d %s -batch "%s"', matlab_executable_path, matlab_command) ;
                % Matlab 2019a-2021a all seem to leak memory when you call getframe() without an
                % X11 server attached.  This is useful as a workaround for that bug.
        else
            bash_command = sprintf('%s -batch "%s"', matlab_executable_path, matlab_command) ;
        end            
        bsub_command = ...
            sprintf('bsub -n %d -eo %s -oo %s %s %s', slot_count, stdouterr_file_name, stdouterr_file_name, options, bash_command) ;
        %fprintf('%s\n', bsub_command) ;
        if isempty(submit_host_name),
            final_command_line = bsub_command ;
        else
            escaped_bsub_command = escape_string_for_bash(bsub_command) ;
            final_command_line = sprintf('ssh -o BatchMode=yes -o StrictHostKeyChecking=no -o ConnectTimeout=20 -q %s %s', ...
                                         submit_host_name, ...
                                         escaped_bsub_command) ;
        end        
        raw_stdout = system_with_error_handling(final_command_line) ;
        stdout = strtrim(raw_stdout) ;  % There are leading newlines and other nonsense in the raw version
        raw_tokens = strsplit(stdout) ;
        is_token_nonempty = cellfun(@(str)(~isempty(str)), raw_tokens) ;
        tokens = raw_tokens(is_token_nonempty) ;
        is_job_token = strcmp(tokens, 'Job') ;
        job_token_index = find(is_job_token, 1) ;
        if isempty(job_token_index) ,
            error('There was a problem submitting the bsub command %s.  Unable to parse output to get job id.  Output was: %s', final_command_line, stdout) ;
        end
        if length(tokens) < job_token_index+3 ,
            error('There was a problem submitting the bsub command %s.  Unable to parse output to get job id.  Output was: %s', final_command_line, stdout) ;
        end
        if ~isequal(tokens{job_token_index+2}, 'is') ,
            error('There was a problem submitting the bsub command %s.  Unable to parse output to get job id.  Output was: %s', final_command_line, stdout) ;
        end
        if ~isequal(tokens{job_token_index+3}, 'submitted') ,
            error('There was a problem submitting the bsub command %s.  Unable to parse output to get job id.  Output was: %s', final_command_line, stdout) ;
        end
        job_id_token = tokens{job_token_index+1} ;
        if length(job_id_token)<2 ,
            error('There was a problem submitting the bsub command %s.  Unable to parse output to get job id.  Output was: %s', final_command_line, stdout) ;
        end
        if ~isequal(job_id_token(1), '<') || ~isequal(job_id_token(end), '>') ,
            error('There was a problem submitting the bsub command %s.  Unable to parse output to get job id.  Output was: %s', final_command_line, stdout) ;
        end        
        job_id_as_string = job_id_token(2:end-1) ;
        job_id = str2double(job_id_as_string) ;
        if ~isfinite(job_id) ,
            error('There was a problem submitting the bsub command %s.  Unable to parse output to get job id.  Output was: %s', final_command_line, stdout) ;
        end
    else
        % Just call the function locally, but use a try/catch to make it more robust.
        if do_try ,
            try
                feval(function_handle, varargin{:}) ;
                job_id = -1 ;  % represents a job that was run locally and exited cleanly
            catch me
                fprintf('Encountered an error while running a local bsub job.  Here''s some information about the error:\n') ;
                fprintf('%s\n', me.getReport()) ;
                job_id = -2 ;  % represents a job that was run locally and errored
            end
        else
            % Run without the try/catch wrapper, usually for debugging
            feval(function_handle, varargin{:}) ;
            job_id = -1 ;  % represents a job that was run locally and exited cleanly
        end
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
    result = '' ;  % fall-through in case of zero args
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
