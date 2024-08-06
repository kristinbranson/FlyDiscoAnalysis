function print_matlab_memory_usage()
    cmd_line = sprintf('pmap %d | tail -n 1 | awk ''{print $2}''', feature('getpid')) ;
    stdout = system_with_error_handling(cmd_line) ;
    fprintf('%s', stdout) ;
end
