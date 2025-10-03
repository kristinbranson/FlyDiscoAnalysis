function result = determineBowlMarkerTemplateFilePath(bowlMarkerTypeThing, raw_plate_id, analysis_protocol_folder_path)
% Determine the absolute file path of the bowl marker template, taking into
% account the raw_plate_id.
% Returns either a full path to a bowl marker template file, or an empty
% array if the given bowlMarkerTypeThing is 'none'.
% This is a pure function.

if ischar(bowlMarkerTypeThing),
  if strcmp(bowlMarkerTypeThing, 'none')
    result = [] ;
  else
    result = fullfile(analysis_protocol_folder_path, bowlMarkerTypeThing) ;
  end
else
  % plate -> bowlmarkertype
  plateids = bowlMarkerTypeThing(1:2:end-1);
  bowlMarkerTypes = bowlMarkerTypeThing(2:2:end);
  if isnumeric(raw_plate_id),
    plateid = num2str(raw_plate_id);
  else
    plateid = raw_plate_id;
  end
  if iscell(plateids)
    i = find(strcmp(num2str(plateid), plateids));
  else
    i = find(str2double(plateid) == plateids,1);
  end
  if isempty(i),
    error('bowlMarkerType not set for plate %s',plateid);
  end
  if strcmp(bowlMarkerTypes{i}, 'none')
    result = [] ;
  else
    result = fullfile(analysis_protocol_folder_path, bowlMarkerTypes{i});
  end
end

end  % function
