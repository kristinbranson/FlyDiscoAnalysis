function result = relpath(path, start_path)
    if ~exist('start_path', 'var') || isempty(start_path) ,
        start_path = pwd() ;
    end
    path_as_object = path_object(path) ;
    start_path_as_object = path_object(start_path) ;
    result_as_path_object = relpath(path_as_object, start_path_as_object) ;
    result = to_char(result_as_path_object) ;
end
