function shadow_list_from_shadowed_function_index = find_shadowed_files(name)
    % name the name of a function whose dependencies you want to search for shadowed
    % functions.
    % Returns a list of all the files that shadow at least one other file of the same name.
    absolute_path_from_function_index = (matlab.codetools.requiredFilesAndProducts(name))' ;
    is_shadowed_from_function_index = cellfun(@is_shadowed_from_absolute_path, absolute_path_from_function_index) ;
    absolute_path_from_shadowed_function_index = absolute_path_from_function_index(is_shadowed_from_function_index) ;
    shadow_list_from_shadowed_function_index = ...
      cellfun(@shadow_list_from_absolute_path, absolute_path_from_shadowed_function_index, 'UniformOutput', false) ;
    % Print them in a nice way
    cellfun(@print_shadow_list, shadow_list_from_shadowed_function_index) ;    
end



function result = is_shadowed_from_absolute_path(function_absolute_path)
    % Whether a function file shadows at least one other function file
    absolute_path_from_instance_index = shadow_list_from_absolute_path(function_absolute_path) ;
    instance_count = length(absolute_path_from_instance_index) ;
    result = (instance_count>1) ;    
end



function absolute_path_from_instance_index = shadow_list_from_absolute_path(function_absolute_path)
    % The list of files that all have a common name, and are all on the path.
    % The first one is the 'winner', i.e. the one that shadows the others/
    [~,function_file_name] = fileparts2(function_absolute_path) ;
    [absolute_path_from_raw_instance_index, info_from_raw_instance_index] = which('-all',function_file_name) ;
    is_relevant_from_raw_instance_index = cellfun(@(str)(isempty(str)||strcmp(str, 'Shadowed')), info_from_raw_instance_index) ;
    absolute_path_from_instance_index = (absolute_path_from_raw_instance_index(is_relevant_from_raw_instance_index)) ;    
end



function print_shadow_list(absolute_path_from_instance_index)
    cellfun(@(absolute_path)(fprintf('%s\n', absolute_path)), absolute_path_from_instance_index) ;
    fprintf('\n') ;
end
