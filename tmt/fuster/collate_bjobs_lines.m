function [numeric_job_status_from_job_index, lsf_status_string_from_job_index] = collate_bjobs_lines(line_from_bjobs_line_index, job_id_from_job_index)
    % If any of the job ids are invalid (maybe b/c they're too old), bjobs will
    % output a line like
    %
    %     Job <42> is not found
    %
    % for that job.  But these all come at the end of the lines for all the valid
    % jobs ids.  (Actually, the good ones go to stdout, the bad to stderr, and
    % Matlab returns them in stdout-then-stderr order in the 2nd output from
    % system().)  This function re-maps the lines in the bjobs output to the job
    % indices.
    line_count = length(line_from_bjobs_line_index) ;
    job_count = length(job_id_from_job_index) ;
    if line_count ~= job_count ,
        error('bjobs line count(%d) is not equal to job_count (%d)', line_count, job_count) ;
    end
    invalid_job_id_from_line_index = nan(size(line_from_bjobs_line_index)) ;
    numeric_job_status_from_line_index = nan(size(line_from_bjobs_line_index)) ;
    lsf_status_string_from_line_index = repmat({''}, size(line_from_bjobs_line_index)) ;
    for line_index = 1 : line_count ,
        bjobs_line = line_from_bjobs_line_index{line_index} ;
            % bjobs_line should be a string like 'DONE', 'EXIT', 'RUN', 'PEND', etc.
            % But it might also be a line like 'Job <42> is not found'.
        [numeric_job_status, lsf_status_string, job_id] = numeric_job_status_from_bjobs_line(bjobs_line) ;
        numeric_job_status_from_line_index(line_index) = numeric_job_status ;
        lsf_status_string_from_line_index{line_index} = lsf_status_string ;
        invalid_job_id_from_line_index(line_index) = job_id ;  % will be nan if this line corresponds to a valid job
    end
    is_valid_from_line_index = isnan(invalid_job_id_from_line_index) ;
    is_invalid_from_line_index = ~is_valid_from_line_index ;
    job_id_from_invalid_line_index = invalid_job_id_from_line_index(is_invalid_from_line_index) ;
    is_invalid_from_job_index = ismember(job_id_from_job_index, job_id_from_invalid_line_index) ;
    is_valid_from_job_index = ~is_invalid_from_job_index ;
    job_id_from_valid_line_index = job_id_from_job_index(is_valid_from_job_index) ;
    job_id_from_line_index = invalid_job_id_from_line_index ;
    job_id_from_line_index(is_valid_from_line_index) = job_id_from_valid_line_index ;  % fill in the job ids for the valid line indices
    line_index_from_job_index = zeros(size(job_id_from_job_index)) ;
    for job_index = 1 : job_count ,
        job_id = job_id_from_job_index(job_index) ;
        line_index = find(job_id_from_line_index==job_id) ;
        if isempty(line_index) ,
            error('Job id %d not found in bjobs output', job_id) ;
        elseif ~isscalar(job_index) ,
            error('Job id %d found multiple times in bjobs output', job_id) ;
        end
        line_index_from_job_index(job_index) = line_index ;
    end        
    numeric_job_status_from_job_index = numeric_job_status_from_line_index(line_index_from_job_index) ;
    lsf_status_string_from_job_index = lsf_status_string_from_line_index(line_index_from_job_index) ;
end
