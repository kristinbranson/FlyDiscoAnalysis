function result = get_single_bsub_job_status(job_id)
    % Possible results are {-1,0,+1}.
    %   -1 means errored out
    %    0 mean running or pending
    %   +1 means completed successfully
    
    if ~isfinite(job_id) ,
        % This means the job has not been submitted yet
        result = nan ;
    elseif job_id == -1 ,
        % This is a job that was run locally and exited cleanly
        result = +1 ;
    elseif job_id == -2 ,
        % This is a job that was run locally and errored
        result = -1 ;
    else
        command_line = sprintf('bjobs %d', job_id) ;
        [status, stdout] = system(command_line) ;
        if status ~= 0 ,
            error('There was a problem running the command %s.  The return code was %d', command_line, status) ;
        end
        lines = strsplit(stdout, '\n') ;
        if length(lines)<2 ,
            error('There was a problem submitting the bsub command %s.  Unable to parse output.  Output was: %s', bsub_command, stdout) ;
        end        
        bjobs_line = lines{2} ;
        tokens = strsplit(bjobs_line) ;
        if length(tokens)<3 ,
            error('There was a problem submitting the bsub command %s.  Unable to parse output.  Output was: %s', bsub_command, stdout) ;
        end
        lsf_status = tokens{3} ;  % Should be string like 'DONE', 'EXIT', 'RUN', 'PEND', etc.    
        if isequal(lsf_status, 'DONE') ,
            result = +1 ;
        elseif isequal(lsf_status, 'EXIT') ,
            % This seems to indicate an exit with something other than a 0 return code
            result = -1 ;
        elseif isequal(lsf_status, 'PEND') || isequal(lsf_status, 'RUN')  || isequal(lsf_status, 'UNKWN'),
            result = 0 ;
        else
            error('Unknown bjobs status string: %s', lsf_status) ;
        end
    end
end
