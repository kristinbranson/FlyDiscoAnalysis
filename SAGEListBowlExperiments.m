function [datamerge,experiment_ids,iswarning] = SAGEListBowlExperiments(varargin)

[dataset,removemissingdata,leftovers] = myparse_nocheck(varargin,...
  'dataset','metadata','removemissingdata',false);

[datamerge,experiment_ids,iswarning] = SAGEGetBowlData(leftovers{:},...
  'dataset',dataset,'removemissingdata',removemissingdata);

% for backwards compatibility
for i = 1:numel(datamerge),
  datamerge(i).line = datamerge(i).line_name;
end
