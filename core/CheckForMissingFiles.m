function [ismissingfile,missingfiles] = CheckForMissingFiles(expdir,dataloc_params,requiredfiles)

% Determines if any of the 'file indicators' listed in requiredfiles is missing from the
% experiment folder.  Note that each 'file indicator' can either be a 
% actual file name or a field name in the dataloc_params struct, the value of
% which contains a file name.  In fact, each file indicator is first checked to
% see whether it is a field in dataloc_params, and if so that field's value is
% used as the file name.  If there is no such field in dataloc_params, the file indicator 
% is treated as a file name.
%
% After each file indicator is resolved to a file name, the existence of the
% file is checked.  On exit, ismissingfile will be true iff some file is missing.
% missingfiles will contain the names of any missing files.

ismissingfile = false;
missingfiles = {};
for i = 1:numel(requiredfiles),
  file_indicator = requiredfiles{i};
  if isfield(dataloc_params,file_indicator),
    fn = dataloc_params.(file_indicator);
  else
    fn = file_indicator ;
  end
  path = fullfile(expdir,fn);
  if ~exist(path,'file'),
    missingfiles{end+1} = path; %#ok<AGROW>
    ismissingfile = true;
  end
end
