function result = read_file_into_cellstring_keeping_newlines(file_name)
    % Read a while file into a cellstring, one line per element.
    % Newlines are discarded.
    fid = fopen(file_name, 'rt') ;
    if fid<0 ,
        error('read_file_into_cell_string_keeping_newlines:unable_to_open_file', 'Unable to open file %s for reading', file_name) ;
    end
    cleaner = onCleanup(@()(fclose(fid))) ;
    result = cell(0,1) ;
    line = fgets(fid) ;
    while ischar(line) ,
        result = vertcat(result, ...
                         {line}) ;  %#ok<AGROW>
        line = fgets(fid) ;
    end
end
