function result = process_rule(rule, name_from_overall_file_index)
    % Take a rule and add the fields overall_file_index_from_input_file_index and overall_file_index_from_output_file_index
    % to it.
    result = rule ;
    name_from_input_file_index = rule.name_from_input_file_index ;
    name_from_output_file_index = rule.name_from_output_file_index ;
    
    function result = overall_file_index_from_file_name(file_name)
        result = find(strcmp(file_name, name_from_overall_file_index)) ;
    end
    
    overall_file_index_from_input_file_index = cellfun(@overall_file_index_from_file_name, name_from_input_file_index) ;
    overall_file_index_from_output_file_index = cellfun(@overall_file_index_from_file_name, name_from_output_file_index) ;
    
    result.overall_file_index_from_input_file_index = overall_file_index_from_input_file_index ;
    result.overall_file_index_from_output_file_index = overall_file_index_from_output_file_index ;    
    
    %result = rmfield(result, 'name_from_input_file_index') ;
    %result = rmfield(result, 'name_from_output_file_index') ;
end
