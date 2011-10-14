function [datamerge,experiment_ids,iswarning] = SAGEListBowlExperiments(varargin)

[data_type,dataset,leftovers] = myparse_nocheck(varargin,...
  'data_type','QuickStats_BackSubStats_meanNConnComps',...
    'dataset','score');

[datamerge,experiment_ids,iswarning] = SAGEGetBowlData(leftovers{:},...
  'dataset',dataset,'data_type',data_type);

% for backwards compatibility
for i = 1:numel(datamerge),
  datamerge(i).line = datamerge(i).line_name;
end
