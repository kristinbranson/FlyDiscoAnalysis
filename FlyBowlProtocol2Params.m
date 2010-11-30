function params = FlyBowlProtocol2Params(paramType,protocolID,varargin)

[protocolFile] = myparse(varargin,...
  'protocolFile','FlyBowlProtocol.txt');
commentChar = '#';

if ~exist(protocolFile,'file'),
  error('Protocol file %s does not exist.',protocolFile);
end
fid = fopen(protocolFile,'r');
if fid < 0,
  error('Could not open protocol file %s for reading.',protocolFile);
end

if ~ischar(protocolID),
  protocolID = num2str(protocolID);
end

didFind = false;
params = {};

% looking for start of parameter type
mode = 0;
while true,
  line = fgetl(fid);
  if ~ischar(line) && line == -1,
    break;
  end
  line = strtrim(line);
  if isempty(line) || line(1) == commentChar,
    continue;
  end
  if mode == 0,
    matches = regexpi(line,'^start (?<paramType>.+),(?<protocolID>[^,]+)$','names','once');
    if isempty(matches),
      error('Protocol parameter type and identifier should be formatted as "start <paramType>,<protocolID>". Instead, got: %s',line);
    end
    paramTypeCurr = matches.paramType;
    protocolIDCurr = matches.protocolID;
    if strcmp(paramTypeCurr,paramType) && strcmp(protocolIDCurr,protocolID),
      didFind = true;
      mode = 1;
    else
      mode = 2;
    end
  else
    if strcmpi(line,['end ',paramTypeCurr]),
      mode = 0;
    end
    if mode == 1,
      matches = regexp(line,'^(?<paramName>[^,]+),(?<paramValue>.+)$','names','once');
      if isempty(matches),
        error('Protocol parameter of type %s could not be parsed: %s',paramTypeCurr,line);
      end
      params(1:2,end+1) = {matches.paramName,eval(matches.paramValue)}; %#ok<AGROW>
    end
  end
end

if ~didFind,
  error('Did not find paramType %s, protocolID %s',paramType,protocolID);
end
