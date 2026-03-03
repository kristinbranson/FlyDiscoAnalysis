function [info,success] = ParseFlyDiscoLog(logfile)

info = struct;
success = false;

fid = fopen(logfile,'r');
if fid < 0,
  error('Could not open file %s',logfile);
end

pattern = 'Running FlyDiscoPipeline.* on (\d{4}-\d{2}-\d{2}) at (\d{2}:\d{2}:\d{2})';
while true,
  s = fgetl(fid);
  if ~ischar(s),
    warning('Could not parse logfile');
    fclose(fid);
    return;
  end
  m = regexp(s,pattern,'once','tokens');
  if ~isempty(m),
    date = m{1};
    time = m{2};
    info.analysis_starttime = [strrep(date,'-',''),'T',strrep(time,':','')];
    break;
  end
end

patterns = {'Matlab version','matlab_version'
  'Source repo','source_repo'
  'Commit hash','commit_hash'};

for i = 1:size(patterns,1),
  while true,
    s = fgetl(fid);
    if ~ischar(s),
      warning('Could not parse logfile');
      fclose(fid);
      return;
    end
    m = regexp(s,patterns{i,1},'once');
    if ~isempty(m),
      s = fgetl(fid);
      if ~ischar(s),
        warning('Could not parse logfile');
        fclose(fid);
        return;
      end
      info.(patterns{i,2}) = strtrim(s);
      break;
    end
  end
end

pattern1 = 'Git status';
pattern2 = 'Git log:';
while true,
  s = fgetl(fid);
  if ~ischar(s),
    warning('Could not parse logfile');
    fclose(fid);
    return;
  end
  m = regexp(s,pattern1,'once');
  if ~isempty(m),
    info.git_status = {};
    while true
      s = fgetl(fid);
      if ~ischar(s),
        warning('Could not parse logfile');
        fclose(fid);
        return;
      end
      if isempty(s),
        continue;
      end
      m = regexp(s,pattern2,'once');
      if ~isempty(m),
        break;
      end
      info.git_status{end+1} = strtrim(s);
    end
    break;
  end
end

info.git_log = {};
while true
  s = fgetl(fid);
  if ~ischar(s),
    warning('Could not parse logfile');
    fclose(fid);
    return;
  end
  m = regexp(s,'^\*','once');
  if isempty(m),
    break;
  end
  info.git_log{end+1} = strtrim(s);
end

pattern = 'Settings values in FlyDiscoPipeline';
while true,
  s = fgetl(fid);
  if ~ischar(s),
    warning('Could not parse logfile');
    fclose(fid);
    return;
  end
  m = regexp(s,pattern,'once');
  if ~isempty(m),
    break;
  end
end

while true,
  s = fgetl(fid);
  if ~ischar(s),
    warning('Could not parse logfile');
    fclose(fid);
    return;
  end
  s = strtrim(s);
  m = regexp(s,'^(.*):(.*)$','once','tokens');
  if isempty(m),
    break;
  end
  m = strtrim(m);
  info.(m{1}) = eval(m{2});
end

pattern = 'Canonical path to analysis protocol folder is';
while true,
  s = fgetl(fid);
  if ~ischar(s),
    warning('Could not parse logfile');
    fclose(fid);
    return;
  end
  m = regexp(s,pattern,'once');
  if ~isempty(m),
    s = fgetl(fid);
    if ~ischar(s),
      warning('Could not parse logfile');
      fclose(fid);
      return;
    end
    info.analysis_protocol_folder = strtrim(s);
  end
  break;
end

fclose(fid);
success = true;

