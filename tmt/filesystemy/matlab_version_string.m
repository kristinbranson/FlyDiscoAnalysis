function result = matlab_version_string()
    % E.g. for Matlab R2022a, returns '2022a'.
    s = ver('Matlab') ;
    full_release_string = s.Release ;  % e.g. '(R2022a)'
    result = full_release_string(3:end-1) ;    
end
