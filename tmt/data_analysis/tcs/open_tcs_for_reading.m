function fid=open_tcs_for_reading(file_name)

fid=fopen(file_name,'r','ieee-le','windows-1252');
if fid<0
  error(sprintf('unable to open file %s',file_name));
end
