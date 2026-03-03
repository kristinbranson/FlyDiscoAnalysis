function [state_from_overall_file_index, state_from_rule_index] = ...
        compute_initial_node_states(name_from_overall_file_index, processed_rule_from_rule_index, ...
                                    is_input_from_overall_file_index_from_rule_index, ...
                                    is_output_from_overall_file_index_from_rule_index, ...
                                    is_a_pure_input_from_overall_file_index)
    overall_file_count = length(name_from_overall_file_index) ;
    rule_count = length(processed_rule_from_rule_index) ;

    state_from_overall_file_index = zeros(overall_file_count, 1) ;
    state_from_rule_index = zeros(rule_count, 1) ;
    
    % Determine the state of the pure input files
    overall_file_index_from_pure_input_index = find(is_a_pure_input_from_overall_file_index) ;
    name_from_pure_input_index = name_from_overall_file_index(is_a_pure_input_from_overall_file_index) ;
    pure_input_count = length(overall_file_index_from_pure_input_index) ;
    for pure_input_index = 1 : pure_input_count ,
        file_name = name_from_pure_input_index{pure_input_index} ;
        does_file_exist = logical(exist(file_name, 'file')) ;
        file_state = 2*does_file_exist-1 ;
        overall_file_index = overall_file_index_from_pure_input_index(pure_input_index) ;
        state_from_overall_file_index(overall_file_index) = file_state ;
    end
   
    was_visited_from_overall_file_index = false(overall_file_count, 1) ;
    was_visited_from_rule_index = false(rule_count, 1) ;
    is_in_current_set_from_overall_file_index = is_a_pure_input_from_overall_file_index ;
    maximum_iteration_count = rule_count ;
    for iteration_index = 1 : maximum_iteration_count ,
        if ~any(is_in_current_set_from_overall_file_index) ,
            break
        end
        % Compute the rules just downstream from the current file set
        is_input_from_current_set_file_index_from_rule_index = ...
            is_input_from_overall_file_index_from_rule_index(is_in_current_set_from_overall_file_index,:) ;
        is_in_current_set_from_rule_index = any(is_input_from_current_set_file_index_from_rule_index, 1)' & ~was_visited_from_rule_index ;
        rule_index_from_current_set_rule_index = find(is_in_current_set_from_rule_index) ;
        current_set_rule_count = length(rule_index_from_current_set_rule_index) ;
        for current_set_rule_index = 1 : current_set_rule_count ,
            rule_index = rule_index_from_current_set_rule_index(current_set_rule_index) ;
            processed_rule = processed_rule_from_rule_index(rule_index) ;
            rule_state = compute_processed_rule_state(processed_rule, state_from_overall_file_index) ;            
            state_from_rule_index(rule_index) = rule_state ;
            if rule_state ~= 0 ,  % if the rule is settled, update the output files to the same state
                % update the state of this rule, and the output files of it
                is_output_from_overall_file_index = is_output_from_overall_file_index_from_rule_index(:,rule_index) ;
                state_from_overall_file_index(is_output_from_overall_file_index) = rule_state ;
            end
        end
        
        % Prep for next iter
        was_visited_from_overall_file_index = was_visited_from_overall_file_index | is_in_current_set_from_overall_file_index ;
        was_visited_from_rule_index = was_visited_from_rule_index | is_in_current_set_from_rule_index ;
        is_in_current_set_from_overall_file_index = (state_from_overall_file_index~=0) & ~was_visited_from_overall_file_index ;
    end    
    if any(is_in_current_set_from_overall_file_index) ,
        error('Internal error: is_in_current_set_from_overall_file_index has nonzero elements after the maximum number of iterations') ;
    end
end
