function [datamerge,experiment_ids,iswarning] = SAGEListBubbleExperiments(varargin)

[dataset,removemissingdata,leftovers] = myparse_nocheck(varargin,...
  'dataset','metadata','removemissingdata',false);

[datamerge,experiment_ids,iswarning] = SAGEGetBubbleData(leftovers{:},...
  'dataset',dataset,'removemissingdata',removemissingdata);

% for backwards compatibility
for i = 1:numel(datamerge),
  datamerge(i).line = datamerge(i).line_name;
end

for i = 1:numel(datamerge)
    datamerge(i).exp_date = datamerge(i).exp_datetime(1:8);
end