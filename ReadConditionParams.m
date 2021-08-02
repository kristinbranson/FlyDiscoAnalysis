function params = ReadConditionParams(filename)

params.flyconditions = {};
params.frameconditions = {};
params.colors = {};
params.markers = {};

fid = fopen(filename,'r');
if fid < 0,
  error('Could not open parameter file %s for reading',filename);
end

while true,
  s = fgetl(fid);
  
  % end of file
  if ~ischar(s), break; end
  
  % comments
  if isempty(s) || ~isempty(regexp(s,'^\s*$','once')) || ...
        ~isempty(regexp(s,'^\s*#','once')),
    continue;
  end
  ss = regexp(strtrim(s),',','split');
  flycondition = '*';
  framecondition = '*';
  marker = 'o';
  isused = false(1,numel(ss));
  for i = 1:min(numel(ss),2),
    m = regexp(ss{i},'^fly:(.*)$','tokens','once');
    if ~isempty(m),
      flycondition = m{1};
      isused(i) = true;
      break;
    end
  end
  
  for i = 1:min(numel(ss),2),
    if isused(i),
      continue;
    end
    m = regexp(ss{i},'^frame:(.*)$','tokens','once');
    if ~isempty(m),
      framecondition = m{1};
      isused(i) = true;
      break;
    end
  end
  i = find(~isused,1);
  if isempty(i),
    continue;
  end
  color = ss{i};
  i = i + 1;
  if i <= numel(ss),
    marker = ss{i};
  end
  colorvec = str2double(regexp(color,'\s+','split'));
  if numel(colorvec) == 3 && ~any(isnan(colorvec)),
    color = colorvec;
  elseif numel(color) ~= 1 || ~ischar(color),
    warning('Could not parse color %s for condition %s, using default black',color,fn);
    color = 'k';
  end
  params.flyconditions{end+1} = flycondition;
  params.frameconditions{end+1} = framecondition;
  params.colors{end+1} = color;
  params.markers{end+1} = marker;
    
end

fclose(fid);