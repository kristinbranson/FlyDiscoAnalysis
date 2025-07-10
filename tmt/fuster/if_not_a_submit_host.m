function result = if_not_a_submit_host(submit_host_name)
    % If the machine we are running is a submit host, return the empty string.
    % Otherwise, return submit_host_name.
    % This function is useful for code that you want to run whether or not you're
    % on a submit host.
    if isempty(submit_host_name) ,
        % If submit_host_name is empty, there's no need to check, so just normalize
        % the input.
        result = '' ;
    else
        % Check if this machine is a submit host
        [return_code, stdouterr] = system('bjobs') ;  %#ok<ASGLU> 
        if return_code == 0 ,
            % bjobs is on the path, so we assume this is a submit host, to we return the
            % empty string.
            result = '' ;
        else
            % bjobs is *not* on the path, so we assume this is *not* a submit host, and
            % pass the argument without modification.
            result = submit_host_name ;
        end
    end
end

