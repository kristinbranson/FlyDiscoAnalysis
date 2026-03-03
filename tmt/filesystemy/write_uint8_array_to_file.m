function write_uint8_array_to_file(file_name, ar) 
    fid = fopen(file_name, 'wb') ;
    if fid<0 , 
        error('Unable to open file %s for writing', file_name) ;
    end
    cleaner = onCleanup(@()(fclose(fid))) ;
    bytes_written_count = fwrite(fid, ar) ;    
    if bytes_written_count ~= length(ar) ,
        error('Unable to write all bytes to file %s', file_name) ;
    end
end
