function result = read_one_file_name_per_line(file_name) 
  fid = fopen(file_name, 'rt') ;
  if fid < 0 , 
    error('Unable to open file %s', file_name) ;
  end
  cleaner = onCleanup(@()(fclose(fid))) ;
  protoresult = textscan(fid, '%s') ;
  if isempty(protoresult) ,
    error('Unable to read file names from %s', file_name) ;
  end
  result = protoresult{1} ;
end
