function [numeric_job_status, lsf_status_string, job_id] = numeric_job_status_from_bjobs_line(bjobs_line)
    % bjobs_line should be string like 'DONE', 'EXIT', 'RUN', 'PEND', etc.
    % Or it can be something like 'Job <42> is not found'.
    
    % Possible numeric_job_status values are {-1,0,+1,nan}.
    %   -1   means errored out
    %    0   mean running or pending
    %   +1   means completed successfully
    %   nan  Job ID not found
    %
    % The lsf_status_string is usually the same as bjobs_line, 
    % unless e.g. lsf_status_string=='Job <42> is not found', then lsf_status_string
    % will be 'NTFND'.
    %
    % The job_id_if_not_found is the job_id if bjobs_line is of the 'Job <n> is not
    % found' form, extracted from bjobs_line.  Otherwise nan.

    if isequal(bjobs_line, 'DONE') ,
        numeric_job_status = +1 ;
        lsf_status_string = bjobs_line ;
        job_id = nan ;
    elseif isequal(bjobs_line, 'EXIT') ,
        % This seems to indicate an exit with something other than a 0 return code
        numeric_job_status = -1 ;
        lsf_status_string = bjobs_line ;
        job_id = nan ;
    elseif isequal(bjobs_line, 'PEND') || isequal(bjobs_line, 'RUN') || isequal(bjobs_line, 'UNKWN') || ...
            isequal(bjobs_line, 'SSUSP') || isequal(bjobs_line, 'PSUSP') || isequal(bjobs_line, 'USUSP') ,
        numeric_job_status = 0 ;
        lsf_status_string = bjobs_line ;
        job_id = nan ;
    else
        if contains(bjobs_line, 'not found') ,
            numeric_job_status = nan;
            lsf_status_string = 'NTFND' ;
            % Need to extract the job ID
            tokens = strsplit(bjobs_line) ;
            if length(tokens) < 2 ,
                error('Unable to parse bjobs line: "%s"', bjobs_line) ;
            end
            job_id_token = tokens{2} ;
            if length(job_id_token)<2 ,
                error('Unable to parse bjobs line: "%s"', bjobs_line) ;
            end
            if ~isequal(job_id_token(1),'<') ,
                error('Unable to parse bjobs line: "%s"', bjobs_line) ;
            end                
            if ~isequal(job_id_token(end),'>') ,
                error('Unable to parse bjobs line: "%s"', bjobs_line) ;
            end
            job_id_as_string = job_id_token(2:end-1) ;
            if isempty(job_id_as_string) ,
                error('Unable to parse bjobs line: "%s"', bjobs_line) ;
            end
            job_id = str2double(job_id_as_string) ;
            if ~isfinite(job_id) || (round(job_id)~=job_id) || job_id<0 ,
                error('Unable to parse bjobs line: "%s"', bjobs_line) ;
            end
        else
            error('Unable to parse bjobs line: "%s"', bjobs_line) ;
        end
    end
end
