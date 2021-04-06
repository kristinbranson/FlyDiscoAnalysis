function [file_names, file_sizes_in_bytes, is_file, is_dir, is_link, mod_time] = list_remote_dir(source_user, source_host, source_path) 
    escaped_source_path = escape_path_for_bash(source_path) ;
    remote_ls_command_line = horzcat('ls -l -A -U -Q --full-time -- ', escaped_source_path) ;
    command_line = horzcat('ssh ', source_user, '@', source_host, ' ', remote_ls_command_line) ;
    [return_code, stdout] = system(command_line) ;
    if return_code ~= 0 ,
        error('list_remote_dir:failed', 'Unable to list the directory %s as user %s on host %s', source_path, source_user, source_host)
    end
    lines_raw = splitlines(strtrim(stdout)) ;
    lines = lines_raw(2:end) ; % drop 1st line
    line_count = length(lines) ;
    file_names = cell(line_count, 1) ;
    file_sizes_in_bytes = zeros(line_count, 1) ;
    is_file = false(line_count, 1) ;
    is_dir = false(line_count, 1) ;
    is_link = false(line_count, 1) ;
    mod_time = NaT(line_count, 1, 'TimeZone', 'UTC') ;
    for i = 1 : line_count , 
        line = lines{i} ;
        [name, size_in_bytes, is_file_this, is_dir_this, is_link_this, mod_time_this] = extract_name_size_and_type_from_ls_long_line(line) ;
        file_names{i} = name ;
        file_sizes_in_bytes(i) = size_in_bytes ;
        is_file(i) = is_file_this ;
        is_dir(i) = is_dir_this ;
        is_link(i) = is_link_this ;
        mod_time(i) = mod_time_this ;
    end
end
