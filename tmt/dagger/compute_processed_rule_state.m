function result = ...
        compute_processed_rule_state(processed_rule, state_from_overall_file_index)
    
    %name_from_input_file_index = rule.name_from_input_file_index ;
    name_from_output_file_index = processed_rule.name_from_output_file_index ;
    overall_file_index_from_input_file_index = processed_rule.overall_file_index_from_input_file_index ;
    %overall_file_index_from_output_file_index = processed_rule.overall_file_index_from_output_file_index ;   
    
    state_from_input_file_index = state_from_overall_file_index(overall_file_index_from_input_file_index) ;
    if any(state_from_input_file_index == -1) ,
        % If there are any bad input files, the rule state is bad
        result = -1 ;
    elseif any(state_from_input_file_index == 0) ,
        % No rules are bad, but some are unsettled
        % That means the rule is unsettled
        result = 0 ;
    else
        % All the input files are good!
        does_exist_from_output_file_index = ...
            cellfun(@(name)(exist(name, 'file')), name_from_output_file_index) ;
        if all(does_exist_from_output_file_index) ,
            % all input files are good
            % all output files exist
            % need to check dates
            mod_time_from_input_file_index = cellfun(@(name)(modtime(name)), overall_file_index_from_input_file_index) ;
            mod_time_from_output_file_index = cellfun(@(name)(modtime(name)), name_from_output_file_index) ;
            earliest_output_mod_time = min(mod_time_from_output_file_index) ;
            most_recent_input_mod_time = max(mod_time_from_input_file_index) ;
            result = (earliest_output_mod_time > most_recent_input_mod_time) ;  % the rule is either good or unsettled, depending
        else
            % all input files exist
            % one or more output files does not exist
            result = 0 ;  % the rule is unsettled
        end
    end    
end
