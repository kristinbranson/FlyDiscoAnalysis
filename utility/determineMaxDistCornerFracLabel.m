function result = determineMaxDistCornerFracLabel(maxDistCornerFracLabel, metadata, field_name)
% Determine the maximum distance from the corner to the registration mark or LED.
% This might depend on bowl.
% This is a pure function.

if numel(maxDistCornerFracLabel) > 1 ,
  plateids = maxDistCornerFracLabel(1:2:end-1);
  cornerfracs = maxDistCornerFracLabel(2:2:end);
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
    error('%s not set for plate %d', field_name, plateid);
  end
  if iscell(cornerfracs)
    result = str2double(cornerfracs{i});
  else
    result = cornerfracs(i);
  end
else
  result = maxDistCornerFracLabel ;
end

end
