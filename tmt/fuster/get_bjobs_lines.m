function result = get_bjobs_lines(job_ids, submit_host_name)
    % Gets the LSF job status for each job_id.
    % Result is a cell array of strings, each string a job status.
    % Job statuses are things like 'DONE', 'EXIT', 'RUN', 'PEND', etc.   
    % In its current incarnation, a better name might be get_LSF_statuses().
    
    if ~exist('submit_host_name', 'var') || isempty(submit_host_name) ,
        submit_host_name = '' ;
    end    
    job_id_count = length(job_ids) ;
    job_id_count_per_call = 10000 ;
    batch_count = ceil(job_id_count / job_id_count_per_call) ;
    bjobs_lines_from_batch_index = cell(batch_count, 1) ;
    for batch_index = 1 : batch_count ,
        first_job_index = job_id_count_per_call * (batch_index-1) + 1 ;
        last_job_index = min(job_id_count_per_call * batch_index, job_id_count) ;
        job_ids_this_batch = job_ids(first_job_index:last_job_index) ;
        job_ids_as_string = sprintf(' %d', job_ids_this_batch) ;
        bjobs_command_line = ['bjobs -o stat -noheader' job_ids_as_string] ;
        if isempty(submit_host_name),
            command_line = bjobs_command_line ;
        else
            escaped_bjobs_command_line = escape_string_for_bash(bjobs_command_line) ;
            command_line = ...
                sprintf('ssh -o BatchMode=yes -o StrictHostKeyChecking=no -o ConnectTimeout=20 -q %s %s', ...
                        submit_host_name, ...
                        escaped_bjobs_command_line) ;
        end
        bjobs_lines = try_bjobs_command_line_a_few_times(command_line, job_ids_this_batch) ;
        bjobs_lines_from_batch_index{batch_index} = bjobs_lines ;
    end
    if isrow(job_ids) && ~isscalar(job_ids) ,
        % Make the shape of the output match the input
        result = vertcat(bjobs_lines_from_batch_index{:})' ;  % row
    else
        result = vertcat(bjobs_lines_from_batch_index{:}) ;  % col
    end
end



function bjobs_lines = try_bjobs_command_line_a_few_times(command_line, job_ids_this_batch)
    maximum_attempt_count = 3 ;
    delay_between_attempts = 10 ;  % s
    did_succeed = false ;
    for attempt_count = 1 : maximum_attempt_count ,
        [~, stdout] = system(command_line) ;
        % Ignore the return code.
        % As of this writing (2023-09-15), bjobs will return a zero/nonzero error code
        % depending only on whether the *last* job ID in the list is known/unknown.
        % This seems like a bug, but whatta ya gonna do?  So we just ignore the return
        % code.  If something has gone very wrong, the parsing of the stdout should
        % fail, and an error will be thrown.
        bjobs_lines = strtrim(strsplit(strtrim(stdout), '\n'))' ;
            % Want a col vector of lines, with no leading or training whitespace on
            % each line.
        if length(bjobs_lines) == length(job_ids_this_batch) ,
            did_succeed = true ;
            break
        end
        warning('tmt:fuster:cant_parse_bjobs_output', ...
                'On attempt %d of %d, there was a problem submitting the bjobs command %s.  Unable to parse output.  Output was: %s', ...
                attempt_count, maximum_attempt_count, command_line, stdout) ;
        pause(delay_between_attempts) ;
    end
    if ~did_succeed ,
        error('tmt:fuster:unable_to_get_bjobs_information', ...
              'All %d attempts to run this bjobs command failed: %s.', maximum_attempt_count) ;
    end
end
