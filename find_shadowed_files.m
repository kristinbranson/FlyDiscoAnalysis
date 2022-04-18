function absolute_path_from_shadowed_function_index = find_shadowed_files(name)
    absolute_path_from_function_index = matlab.codetools.requiredFilesAndProducts(name) ;
    is_shadowed_from_function_index = cellfun(@is_shadowed_from_absolute_path, absolute_path_from_function_index) ;
    absolute_path_from_shadowed_function_index = absolute_path_from_function_index(is_shadowed_from_function_index) ;
end



function result = is_shadowed_from_absolute_path(function_absolute_path)
    [~,function_file_name] = fileparts2(function_absolute_path) ;
    [absolute_path_from_raw_instance_index, info_from_raw_instance_index] = which('-all',function_file_name) ;
    is_relevant_from_raw_instance_index = cellfun(@(str)(isempty(str)||strcmp(str, 'Shadowed')), info_from_raw_instance_index) ;
    absolute_path_from_instance_index = absolute_path_from_raw_instance_index(is_relevant_from_raw_instance_index) ;    
    instance_count = length(absolute_path_from_instance_index) ;
    if instance_count>1 ,
        absolute_path_from_instance_index
    end
    result = (instance_count>1) ;    
end
