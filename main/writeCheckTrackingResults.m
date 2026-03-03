function writeCheckTrackingResults(savefile,isok,isredo,flag,notes)

fid = fopen(savefile,'w');
if fid < 0,
  error('Could not open file %s for writing CheckTrackingResults',savefile);
end
fprintf(fid,'ok:%d\n',isok);
fprintf(fid,'redo:%d\n',isredo);
fprintf(fid,'flag:%s\n',flag);
if ischar(notes),
  notes = cellstr(notes);
end
fprintf(fid,'notes:%s\n',notes{1});
if length(notes) > 1,
  fprintf(fid,'%s\n',notes{2:end});
end
fclose(fid);
