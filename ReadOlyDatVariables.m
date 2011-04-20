function stats = ReadOlyDatBehaviorVariables(file)

if ~exist(file,'file'),
  error('OlyDatBehaviorVariables file %s does not exist',file);
end
fid = fopen(file,'r');

stats = {};
while true,
  s = fgetl(fid);
  if ~ischar(s), break; end
  s = strtrim(s);
  if isempty(s),
    continue;
  end
  if ~isempty(regexp(s,'^\s*#','once')),
    continue;
  end
  stats{end+1} = s; %#ok<AGROW>
end

fclose(fid);