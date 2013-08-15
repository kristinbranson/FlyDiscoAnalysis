function params = ReadStatsPlotParams(paramsfile)

if ~exist(paramsfile,'file'),
  error('Stats plot params file %s does not exist',paramsfile);
end
fid = fopen(paramsfile,'r');
if fid < 1,
  error('Could not open file %s for reading',paramsfile);
end


params = struct;
current_group = '';
while true,
  
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end
  s = strtrim(s);
  
  if isempty(s) || s(1) == '#'
    continue;
  end
  
  if s(1) == '*',
    current_group = s(2:end);
    if isempty(current_group),
      error('group read as empty: >%s<',s);
    end
    params.(current_group) = {};
    continue;
  end

  if isempty(current_group),
    error('incorrect formatting for file %s, should start with group',paramsfile);
  end

  params.(current_group){end+1} = s;  
  
end

fclose(fid);