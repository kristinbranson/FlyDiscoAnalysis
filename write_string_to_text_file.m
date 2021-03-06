function write_string_to_text_file(file_name, string)
    fid = fopen(file_name, 'wt') ;
    if fid < 0,
        error('Unable to open file %s for writing', file_name) ;
    end
    cleaner = onCleanup(@()(fclose(fid))) ;
    fprintf(fid, '%s', string) ;        
end
