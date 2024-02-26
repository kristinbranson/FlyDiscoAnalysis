function result = determineMaxDistCornerFracLabel(maxDistCornerFracLabel, plate)
% Determine the maximum distance from the corner to the registration mark.
% This might depend on plate.
% This is a pure function.

if numel(maxDistCornerFracLabel) > 1 ,
  plateids = maxDistCornerFracLabel(1:2:end-1);
  cornerfracs = maxDistCornerFracLabel(2:2:end);
  if isnumeric(plate),
    plateid = num2str(plate);
  else
    plateid = plate;
  end
  if iscell(plateids)
    i = find(strcmp(num2str(plateid), plateids));
  else
    i = find(str2double(plateid) == plateids,1);
  end
  if isempty(i),
    error('%s not set for plate %d', 'maxDistCornerFrac_BowlLabel', plateid);
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
