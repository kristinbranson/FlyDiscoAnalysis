function [state_from_rule_index, state_from_overall_file_index] = ...
    dagger(rule_from_rule_index, maximum_running_slot_count, do_actually_submit, maximum_wait_time, do_show_progress_bar)
    % Possible states are {-1,0,+1}.
    %    -1 means bad
    %     0 means unsettled
    %    +1 means good
    if ~exist('maximum_running_slot_count', 'var') || isempty(maximum_running_slot_count) ,
        maximum_running_slot_count = inf ;
    end
    if ~exist('do_actually_submit', 'var') || isempty(do_actually_submit) ,
        do_actually_submit = true ;
    end
    if ~exist('maximum_wait_time', 'var') || isempty(maximum_wait_time) ,
        maximum_wait_time = inf ;
    end
    if ~exist('do_show_progress_bar', 'var') || isempty(do_show_progress_bar) ,
        do_show_progress_bar = true ;
    end

    % Get the indices of the rules, files involved
    [name_from_overall_file_index, ...
     processed_rule_from_rule_index, ...
     is_input_from_overall_file_index_from_rule_index, ...
     is_output_from_overall_file_index_from_rule_index, ...
     is_downstream_from_rule_index_from_rule_index, ...
     is_downstream_from_overall_file_index_from_rule_index, ...
     is_a_pure_input_from_overall_file_index] = ...
        compute_rule_file_lookup_tables(rule_from_rule_index)
    
    % Get the number of files, rules
    %overall_file_count = length(name_from_overall_file_index) ;
    rule_count = length(rule_from_rule_index) ;    
    
    % Compute the initial states of all the files, all the rules
    % Each node (a rule or file) can be good (+1), bad (-1), or unsettled (0).
    % Everything starts out unsettled.
    % Dagger exits when there are no unsettled rules left.
    % Rules become good my running and exiting with no error.
    % Rules become bad by running and exiting with an error, or
    % by being downstream from a bad rule, or by being downstream from a bad file.
    % Files become good by being an output of a good rule, or by being a 'pure' input (not being an output
    % of any rule) and existing.
    % Files become bad by being an output of a bad rule, or by being a pure input
    % and not existing.
    % A node that is either good or bad is 'settled'.

    [state_from_overall_file_index, state_from_rule_index] = ...
        compute_initial_node_states(name_from_overall_file_index, processed_rule_from_rule_index, ...
                                    is_input_from_overall_file_index_from_rule_index, ...
                                    is_output_from_overall_file_index_from_rule_index, ...
                                    is_a_pure_input_from_overall_file_index) ;
                               
   
    
    
    is_time_up = false ;
    is_rule_settled_from_rule_index = (state_from_rule_index~=0) ;
    % is_good_if_an_input_from_overall_file_index_from_rule_index = ...
    %     ~is_input_from_overall_file_index_from_rule_index | (state_from_overall_file_index==1) ;
    % all_inputs_are_good_from_rule_index = all(is_good_if_an_input_from_overall_file_index_from_rule_index, 1)' ;
    % is_runnable_from_rule_index = (state_from_rule_index==0) & all_inputs_are_good_from_rule_index ;  % means rule is ready-to-run
    %     % means rule is unsettled (0) and means inputs are good (+1)
    settled_rule_count = sum(is_rule_settled_from_rule_index) ;
    job_id_from_rule_index = nan(rule_count,1) ;  % nan means not yet submitted
    job_status_from_rule_index = nan(rule_count,1) ;  % nan means the job_id is nan
    was_submitted_from_rule_index = false(rule_count, 1) ;  % this stays true even after the job exits
    is_running_from_rule_index = false(rule_count, 1) ;  % this goes false gain after the job exits
    slot_count_from_rule_index = [rule_from_rule_index.slot_count]' ;
    is_progress_still_possible = true ;
    if do_show_progress_bar ,
        progress_bar = progress_bar_object(rule_count) ;
        progress_bar.update(settled_rule_count) ;
    end
    last_settled_rule_count = settled_rule_count ;
    ticId = tic() ;
    while settled_rule_count<rule_count && ~is_time_up && is_progress_still_possible ,        
        % Update what jobs are running now, and their job statuses
        was_running_from_rule_index = is_running_from_rule_index ;
        job_id_from_was_running_index = job_id_from_rule_index(was_running_from_rule_index) ;
        job_status_from_was_running_index = get_bsub_job_status(job_id_from_was_running_index) ; 
        job_status_from_rule_index(was_running_from_rule_index) = job_status_from_was_running_index ;        
        is_running_from_rule_index = was_submitted_from_rule_index & (job_status_from_rule_index==0) ;
        
        % Update state_from_rule_index, state_from_overall_file_index
        did_job_just_exit_from_rule_index = was_running_from_rule_index & ~is_running_from_rule_index ;
        rule_index_from_just_exited_job_index = find(did_job_just_exit_from_rule_index) ;
        just_exited_job_count = length(rule_index_from_just_exited_job_index) ;
        for just_exited_job_index = 1 : just_exited_job_count ,
            rule_index = rule_index_from_just_exited_job_index(just_exited_job_index) ;
            rule_state = job_status_from_rule_index(rule_index) ;
            if rule_state == +1 ,
                state_from_rule_index(rule_index) = +1 ;  % update the rule state
                % TODO: If the rule_state is -1, shouldn't we update all downstream rules and
                % files?
                is_output_from_overall_file_index = is_output_from_overall_file_index_from_rule_index(:,rule_index) ;
                state_from_overall_file_index(is_output_from_overall_file_index) = +1 ;  % update the file states for the rule outputs
            elseif rule_state == -1 ,
                % Mark all the downstream nodes as bad
                state_from_rule_index(rule_index) = -1 ;  % update the rule state
                is_downstream_from_rule_index = is_downstream_from_rule_index_from_rule_index(:,rule_index) ;
                state_from_rule_index(is_downstream_from_rule_index) = -1 ;
                is_downstream_from_overall_file_index = is_downstream_from_overall_file_index_from_rule_index(:,rule_index) ;
                state_from_overall_file_index(is_downstream_from_overall_file_index) = -1 ;
            else
                error('Internal error: A just-exited job has a job status of %d', rule_state) ;
            end
        end
        
        % Update what rules are settled
        is_rule_settled_from_rule_index = (state_from_rule_index~=0) ;
        settled_rule_count = sum(is_rule_settled_from_rule_index) ;
        newly_settled_rule_count = settled_rule_count - last_settled_rule_count ;
        
        % Update the progress bar
        if do_show_progress_bar ,
            progress_bar.update(newly_settled_rule_count) ;
        end
        
        % Update which rules are runnable
        % A rule is "done" if the job for it exited, or it did not need to be run
        is_good_if_an_input_from_overall_file_index_from_rule_index = ...
            ~is_input_from_overall_file_index_from_rule_index | (state_from_overall_file_index==1) ;
        all_inputs_are_good_from_rule_index = all(is_good_if_an_input_from_overall_file_index_from_rule_index, 1)' ;
        is_runnable_from_rule_index = (state_from_rule_index==0) & all_inputs_are_good_from_rule_index ;
            % means rule is unsettled (0) and means inputs are good (+1)
            % means rule is ready-to-run
        
        % Compute how many slots are in use, and thus how many slots we can use for new
        % job submissions
        carryover_slot_count = sum(slot_count_from_rule_index(is_running_from_rule_index)) ;
        maximum_new_slot_count = maximum_running_slot_count - carryover_slot_count ;
        
        % Submit new jobs, if appropriate
        is_submittable_from_rule_index = ~was_submitted_from_rule_index & is_runnable_from_rule_index ;
            % to be "submittable", a rule has to be runnable and not-yet-submitted
        rule_index_from_submittable_index = find(is_submittable_from_rule_index) ;
        slot_count_from_submittable_index = slot_count_from_rule_index(rule_index_from_submittable_index) ;
        will_submit_from_submittable_index = determine_which_jobs_to_submit(slot_count_from_submittable_index, maximum_new_slot_count) ;
        rule_index_from_submission_index = rule_index_from_submittable_index(will_submit_from_submittable_index) ;
        submission_count = length(rule_index_from_submission_index) ;
        for submission_index = 1 : submission_count ,
            rule_index = rule_index_from_submission_index(submission_index) ;
            rule = rule_from_rule_index(rule_index) ;
            function_handle = rule.function_handle ;
            other_arguments = rule.other_arguments ;
            slot_count = rule.slot_count ;
            stdouterr_file_name = rule.stdouterr_file_name ;
            bsub_option_string = rule.bsub_option_string ;            
            do_use_xvfb = rule.do_use_xvfb ;
            job_id = ...
                bsub(do_actually_submit, ...
                     slot_count, ...
                     stdouterr_file_name, ...
                     bsub_option_string, ...
                     do_use_xvfb, ...
                     function_handle, ...
                     other_arguments{:}) ;
            job_id_from_rule_index(rule_index) = job_id ;
            was_submitted_from_rule_index(rule_index) = true ;
            is_running_from_rule_index(rule_index) = true ;  % for local jobs, this is not exactly true, but it works.
        end        
        
        % If no jobs are in progress at this point, then no future
        % progress can happen
        is_progress_still_possible = any(is_running_from_rule_index) ;

        % See if our time is up
        is_time_up = (toc(ticId) > maximum_wait_time) ;
       
        % Update this
        last_settled_rule_count = settled_rule_count ;
        
        % If we're actually submitting to the cluster, pause here a bit so we
        % don't submit too too many bjobs commands
        if do_actually_submit ,
            pause(1) ;
        end        
    end
%     if do_show_progress_bar ,
%         progress_bar.finish_up() ;
%     end
end
