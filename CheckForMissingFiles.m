function [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles)

ismissingfile = false;
missingfiles = {};
for i = 1:numel(requiredfiles),
  fn = requiredfiles{i};
  if isfield(dataloc_params,fn),
    fn = dataloc_params.(fn);
  end
  fn = fullfile(expdir,fn);
  if ~exist(fn,'file'),
    missingfiles{end+1} = fn; %#ok<AGROW>
    ismissingfile = true;
  end
end