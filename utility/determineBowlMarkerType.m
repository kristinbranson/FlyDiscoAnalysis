function result = determineBowlMarkerType(bowlMarkerType, metadata, settingsdir, analysis_protocol)
% Determine the bowl marker type. 
% This is a pure function.

result = bowlMarkerType ;
if ischar(bowlMarkerType),
  if ~ismember(bowlMarkerType,{'gradient'}),
    result = fullfile(settingsdir,analysis_protocol,bowlMarkerType);
  end
else
  % plate -> bowlmarkertype
  plateids = bowlMarkerType(1:2:end-1);
  bowlmarkertypes = bowlMarkerType(2:2:end);
  if isnumeric(metadata.plate),
    plateid = num2str(metadata.plate);
  else
    plateid = metadata.plate;
  end
  if iscell(plateids)
    i = find(strcmp(num2str(plateid), plateids));
  else
    i = find(str2double(plateid) == plateids,1);
  end
  if isempty(i),
    error('bowlMarkerType not set for plate %s',plateid);
  end
  if ~ismember(bowlmarkertypes{i},{'gradient'}),
    result = fullfile(settingsdir,analysis_protocol,bowlmarkertypes{i});
  end
end

end
