function [name_from_overall_file_index, ...
          processed_rule_from_rule_index, ...
          is_input_from_overall_file_index_from_rule_index, ...
          is_output_from_overall_file_index_from_rule_index, ...
          is_downstream_from_rule_index_from_rule_index, ...
          is_downstream_from_overall_file_index_from_rule_index, ...
          is_a_pure_input_from_overall_file_index] = ...
        compute_rule_file_lookup_tables(rule_from_rule_index)
    % From the rules, computes a bunch of structures that will be useful 
    % when we actually want to execute the DAG defined by the rules.    
    
    % First compile a list of all the unique files in all the rules
    name_from_overall_input_file_index = unique([rule_from_rule_index.name_from_input_file_index]) ;
    name_from_overall_output_file_index = unique([rule_from_rule_index.name_from_output_file_index]) ;
    name_from_overall_file_index = unique(horzcat(name_from_overall_input_file_index, name_from_overall_output_file_index)) ;
    
    % Make a closure that closes over name_from_overall_file_index
    function result = process_rule_wrapper(rule)
        result = process_rule(rule, name_from_overall_file_index) ;
    end
    
    % Process the rules so we have the overall file indices for each rule
    processed_rule_from_rule_index = arrayfun(@process_rule_wrapper, rule_from_rule_index) ;
    
    % Make the arrays that allow quick lookups of file given rule or rule given file
    overall_file_count = length(name_from_overall_file_index) ;
    rule_count = length(rule_from_rule_index) ;
    is_input_from_overall_file_index_from_rule_index = false(overall_file_count, rule_count) ;
    is_output_from_overall_file_index_from_rule_index = false(overall_file_count, rule_count) ;
    for rule_index = 1 : rule_count , 
        processed_rule = processed_rule_from_rule_index(rule_index) ;
        overall_file_index_from_input_file_index = processed_rule.overall_file_index_from_input_file_index ;
        overall_file_index_from_output_file_index = processed_rule.overall_file_index_from_output_file_index ;
        is_input_from_overall_file_index_from_rule_index(overall_file_index_from_input_file_index, rule_index) = true ;        
        is_output_from_overall_file_index_from_rule_index(overall_file_index_from_output_file_index, rule_index) = true ;        
    end
    
    % Compute the "downstream" matrix, allows us to quickly determine if a rule is
    % downstream from another rule.
    RR = is_input_from_overall_file_index_from_rule_index' * is_output_from_overall_file_index_from_rule_index + eye(rule_count) ;
       % rule_count x rule_count, RR(i,j) indicates whether rule i is immediately
       % downstream of rule j, where we count each rule as being immediately downstream of itself. 
    D_putative_old = [] ;
    D_putative = RR ;
    while ~isequal(D_putative, D_putative_old) ,
        D_putative_old = D_putative ;
        D_putative = min(RR*D_putative, ones(rule_count)) ;
    end
    is_downstream_from_rule_index_from_rule_index = D_putative ;  % rule_count x rule_count
    
    % Compute the other downstream matrix
    % F x R
    is_downstream_from_overall_file_index_from_rule_index = ...
        min(is_output_from_overall_file_index_from_rule_index * is_downstream_from_rule_index_from_rule_index, ones(overall_file_count, rule_count)) ;
    
    % Compute which files are "pure" inputs
    is_an_input_from_overall_file_index = max(is_input_from_overall_file_index_from_rule_index, [], 2) ;
    is_an_output_from_overall_file_index = max(is_output_from_overall_file_index_from_rule_index, [], 2) ;
    is_a_pure_input_from_overall_file_index = is_an_input_from_overall_file_index & ~is_an_output_from_overall_file_index ;
end
