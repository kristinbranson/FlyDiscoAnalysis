function [stats,success,errmsg] = readUFMFDiagnostics(filename)

success = false;
errmsg = '';
stats = struct;

% open the file
fid = fopen(filename,'r');
if fid < 0,
  errmsg = sprintf('Could not open file %s for reading',filename);
  return;
end

while true,
  
  l = ftell(fid);
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end
  fseek(fid,l,'bof');
  
  switch s,
    case 'streamStart'
      [stats.stream,success1,errmsg] = readUFMFStream(fid);
      if ~success1,
        fclose(fid);
        return;
      end
    case 'summaryStart'
      [stats.summary,success1,errmsg] = readUFMFSummary(fid);
      if ~success1,
        fclose(fid);
        return;
      end
    otherwise
      errmsg = sprintf('Unknown token "%s".',s);
      fclose(fid);
      return;
  end
  
end

success = true;
fclose(fid);

function [stream,success,errmsg] = readUFMFStream(fid)

success = false;
errmsg = '';
stream = struct;

% read streamStart
s = fgetl(fid);
if ~strcmp(s,'streamStart'),
  errmsg = 'Stream does not start with the keyword "streamStart" on its own line';
  return;
end

% read headers
s = fgetl(fid);
if ~ischar(s),
  errmsg = 'No headers found.';
  return;
end
% split at ,
headers = regexp(s,',','split');
% remove whitespace
for i = 1:length(headers),
  headers{i} = makeLegalName(strtrim(headers{i}));
end
nfields = length(headers);

% loop over frames
line = 1;
while true,
  % read one line of data
  s = fgetl(fid);
  if ~ischar(s),
    errmsg = 'Did not find "streamEnd" token.';
    return;
  end
  if strcmp(s,'streamEnd'),
    break;
  end
  % split at ,
  datas = regexp(s,',','split');
  if length(datas) ~= nfields,
    errmsg = sprintf('Error reading line %d of the stream: number of fields in this line %d does not match number of headers %d',line,length(datas),nfields);
    return;
  end
  % try converting to number
  data = str2double(datas);
  % loop over each field
  for i = 1:length(data),
    fn = headers{i};
    stream.(fn)(line) = data(i);
  end
  line = line + 1;
end

success = true;

function [summary,success,errmsg] = readUFMFSummary(fid)

success = false;
errmsg = '';
summary = struct;

% read summarySummary
s = fgetl(fid);
if ~strcmp(s,'summaryStart'),
  errmsg = 'Summary does not start with the keyword "summaryStart" on its own line';
  return;
end

% read headers
s = fgetl(fid);
if ~ischar(s),
  errmsg = 'No headers found.';
  return;
end
% split at ,
headers = regexp(s,',','split');
% remove whitespace
for i = 1:length(headers),
  headers{i} = makeLegalName(strtrim(headers{i}));
end
nfields = length(headers);

% read one line of data
s = fgetl(fid);
if ~ischar(s),
  errmsg = 'Did not find summary data.';
  return;
end
% split at ,
datas = regexp(s,',','split');
if length(datas) ~= nfields,
  errmsg = sprintf('Error reading summary data line: number of fields in this line %d does not match number of headers %d',length(datas),nfields);
  return;
end
% try converting to number
data = str2double(datas);
% loop over each field
for i = 1:length(data),
  fn = headers{i};
  summary.(fn) = data(i);
end

s = fgetl(fid);
if ~ischar(s) || ~strcmp(s,'summaryEnd'),
  errmsg = 'Did not find summaryEnd token.';
  return;
end

success = true;

function s = makeLegalName(s)

s = regexprep(s,'<=','_leq_');
s = regexprep(s,'>=','_geq_');
s = regexprep(s,'<','_lt_');
s = regexprep(s,'>','_gt_');
s = regexprep(s,'\.','_');

s = regexprep(s,'[^A-Za-z_0-9]','');
