function [isok,isredo,flag,notes,didread] = readCheckTrackingResults(savefile)

notes = '';
isok = false;
isredo = false;
flag = 'None';
didread = false;
if ~exist(savefile,'file'),
  return;
end

fid = fopen(savefile,'r');
if fid < 0,
  return;
end

while true,
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end
  names = regexp(s,'^(?<var>[^:]+):(?<val>.*)$','names','once');
  if isempty(names),
    continue;
  end
  switch names.var,
    case 'ok',
      isok = str2double(names.val)~=0;
    case 'redo',
      isredo = str2double(names.val)~=0;
    case 'flag',
      flag = names.val;
    otherwise,
      notes = names.val;
      isfirst = true;
      while true,
        s = fgetl(fid);
        if ~ischar(s), break; end
        if isfirst,
          notes = {notes,s};
          isfirst = false;
        else
          notes = [notes,{s}];
        end
      end
      break;
  end
end
fclose(fid);

didread = true;