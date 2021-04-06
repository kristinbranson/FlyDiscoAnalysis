function [name, size_in_bytes, is_file, is_dir, is_link, mod_time] = extract_name_size_and_type_from_ls_long_line(line)
    % We assume line looks like this: '-rw-r--r--  1 taylora scicompsoft     278 2020-12-02 17:09:49.027303272 -0500 "test_bw_smooth.m"'
    tokens = strsplit(line) ;
    size_in_bytes = str2double(tokens{5}) ;
    parts = strsplit(line, '"') ;
    name = parts{2} ;
    file_type_char = line(1) ;
    is_file = isequal(file_type_char, '-') ;
    is_dir =  isequal(file_type_char, 'd') ;
    is_link = isequal(file_type_char, 'l') ;
    mod_date_as_string = tokens{6} ;
    mod_time_as_string_with_ns = tokens{7} ;
    mod_time_as_string = mod_time_as_string_with_ns(1:15) ;  % only out to ms
    utc_offset_as_string = tokens{8} ;
    mod_time_as_string = horzcat(mod_date_as_string, ' ', mod_time_as_string, ' ', utc_offset_as_string) ;  % date, time
    time_format = 'yyyy-MM-dd HH:mm:ss.SSSSSS XX' ;
    mod_time = datetime(mod_time_as_string, 'InputFormat', time_format, 'TimeZone', 'UTC') ;  % Want a datetime in UTC timezone 
end
