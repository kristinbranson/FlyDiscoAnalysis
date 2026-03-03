function touch(file_name) 
    fid = fopen(file_name, 'a') ;
    if fid<0 , 
        error('Unable to open file %s for touching', file_name) ;
    end
    fclose(fid) ;
end
