function result = determineMaxDistCornerFracBowlLabel(maxDistCornerFrac_BowlLabel, metadata)
% Determine the maximum distance from the corner to the registration mark.
% This might depend on bowl.
% This is a pure function.

if numel(maxDistCornerFrac_BowlLabel) > 1 ,
  plateids = maxDistCornerFrac_BowlLabel(1:2:end-1);
  cornerfracs = maxDistCornerFrac_BowlLabel(2:2:end);
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
    error('maxDistCornerFrac_BowlLabel not set for plate %d',plateid);
  end
  if iscell(cornerfracs)
    result = str2double(cornerfracs{i});
  else
    result = cornerfracs(i);
  end
else
  result = maxDistCornerFrac_BowlLabel ;
end

end
