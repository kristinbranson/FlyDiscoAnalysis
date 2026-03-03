function [numeric_job_status_from_job_index, lsf_status_string_from_job_index] = ...
        get_bsub_job_status(job_id_from_job_index, submit_host_name)
    % Get the status of each job_ids
    % Possible results are {-1,0,+1,nan}.
    %   -1   means errored out
    %    0   means running or pending
    %   +1   means completed successfully
    %   nan  means job ID not found
    %
    % The seconds return value gives the raw LSF status for each job, as a
    % string.  LSF statuses are things like 'DONE', 'EXIT', 'RUN', 'PEND', etc. 

    if ~exist('submit_host_name', 'var') || isempty(submit_host_name) ,
        submit_host_name = '' ;
    end    
    numeric_job_status_from_job_index = nan(size(job_id_from_job_index)) ;
    lsf_status_string_from_job_index = repmat({''}, size(job_id_from_job_index)) ;
    is_not_yet_submitted_from_job_index = isnan(job_id_from_job_index) ;
    was_run_locally_from_job_index = (job_id_from_job_index<0) ;  % means the job was run locally
    was_run_locally_and_exited_cleanly_from_job_index = (job_id_from_job_index==-1) ;
    was_run_locally_and_errored_from_job_index = (job_id_from_job_index==-2) ;
    numeric_job_status_from_job_index(was_run_locally_and_exited_cleanly_from_job_index) = +1 ;
    numeric_job_status_from_job_index(was_run_locally_and_errored_from_job_index) = -1 ;
    if all(is_not_yet_submitted_from_job_index | was_run_locally_from_job_index) ,
        return
    end
    was_submitted_from_job_index = ~(was_run_locally_from_job_index | is_not_yet_submitted_from_job_index) ;
    job_id_from_submitted_job_index = job_id_from_job_index(was_submitted_from_job_index) ;
    % We don't support repeated job ids (unless the jobs are local), so check for that
    if length(unique(job_id_from_submitted_job_index)) ~= length(job_id_from_submitted_job_index) ,
        error('tmt:fuster:repeated_job_ids_not_supported', 'Repeated job ids not supported') ;
    end
    line_from_bjob_line_index = get_bjobs_lines(job_id_from_submitted_job_index, submit_host_name) ;
    [numeric_job_status_from_submitted_job_index, lsf_status_string_from_submitted_job_index] = ...
        collate_bjobs_lines(line_from_bjob_line_index, job_id_from_submitted_job_index) ;
    numeric_job_status_from_job_index(was_submitted_from_job_index) = numeric_job_status_from_submitted_job_index ;
    lsf_status_string_from_job_index(was_submitted_from_job_index) = lsf_status_string_from_submitted_job_index ;
end
