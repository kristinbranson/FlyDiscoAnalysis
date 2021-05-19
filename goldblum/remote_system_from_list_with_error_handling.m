function stdout = remote_system_from_list_with_error_handling(user_name, host_name, remote_command_line_as_list)
    % Run the system command, but taking a list of tokens rather than a string, and
    % running on a remote host.  Uses ssh, which needs to be set up for passowrdless
    % login as the indicated user.
    % Each element of command_line_as_list is escaped for bash, then composed into a
    % single string, then submitted to system_with_error_handling().
    
    % Escape all the elements of command_line_as_list
    escaped_remote_command_line_as_list = cellfun(@escape_string_for_bash, remote_command_line_as_list, 'UniformOutput', false) ;
    
    % Build up the command line by adding space between elements
    remote_command_line = space_out(escaped_remote_command_line_as_list) ;

    % Command line
    command_line_as_list = {'ssh', '-l', user_name, host_name, remote_command_line} ; 
    
    % Actually run the command
    stdout = system_from_list_with_error_handling(command_line_as_list) ;
end



function result = space_out(list)
    result = '' ;
    count = length(list) ;
    for i = 1 : count ,
        if i==1 ,
            result = list{i} ;
        else
            result = [result ' ' list{i}] ;  %#ok<AGROW>
        end 
    end
end