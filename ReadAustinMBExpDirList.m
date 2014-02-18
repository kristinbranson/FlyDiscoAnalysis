function expdirs = ReadAustinMBExpDirList(filename)

 expdirs = {};
 fid = fopen(filename,'r');
 while true,
   s = fgetl(fid);
   if ~ischar(s),
     break;
   end
   s = strtrim(s);
   if isempty(s),
     continue;
   end
   s = strrep(s,'/Volumes','/groups/branson');
   expdirs{end+1} = s;
 end
 fclose(fid);