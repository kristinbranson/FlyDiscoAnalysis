function result = read_file_into_uint8_array(file_name)
    fid = fopen(file_name, 'rb') ;
    if fid<0 ,
        error('Unable to open file %s for reading', file_name) ;
    end
    cleaner = onCleanup(@()(fclose(fid))) ;
    result = fread(fid, inf, 'uint8=>uint8') ;
end
