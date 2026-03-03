function result = compute_does_shadow_a_builtin(m_file_path) 
    [~,function_name] = fileparts(m_file_path) ;
    [file_path_from_index, details_from_index] = which(function_name, '-all') ;
    is_shadowed_from_index = strcmpi(details_from_index, 'shadowed') ;
    is_self_from_index = strcmpi(file_path_from_index, m_file_path) ;
    is_shadowed_and_not_self_from_index = is_shadowed_from_index&~is_self_from_index ;
    if any(is_shadowed_and_not_self_from_index) ,
        % file_path_from_index
        % details_from_index
        file_path_from_shadowed_index = file_path_from_index(is_shadowed_from_index) ;
        fprintf('Library file %s shadows:\n', m_file_path) ;
        for i = 1 : length(file_path_from_shadowed_index) ,
            fprintf('    %s\n', file_path_from_shadowed_index{i}) ;
        end
        fprintf('\n') ;        
        result = true ;
    else
        result = false ;
    end
end
