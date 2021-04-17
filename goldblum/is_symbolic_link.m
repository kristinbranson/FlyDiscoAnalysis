function result = is_symbolic_link(path)
    escaped_path = escape_string_for_bash(path) ;
    command_line = sprintf('test -L %s', escaped_path) ;
    retval = system(command_line) ;
    result = (retval==0) ;
end
