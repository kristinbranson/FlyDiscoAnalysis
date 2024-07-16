function [ismissingfile, missingfiles] = CheckForMissingFilesExplicit(expdir, requiredfiles)

% Determines if any of the file paths listed in requiredfiles is missing from the
% experiment folder.  Each element of requiredfiles is assumed to be a simple
% file path.  (Constrast with CheckForMissingFiles().)
%
% The existence of the file is checked.  On exit, ismissingfile will be true
% iff some file is missing. missingfiles will contain the names of any missing
% files.

ismissingfile = false ;
missingfiles = cell(1,0) ;
for i = 1:numel(requiredfiles) ,
  fn = requiredfiles{i} ;
  path = fullfile(expdir,fn) ;
  if ~exist(path, 'file') ,
    missingfiles{end+1} = path ; %#ok<AGROW>
    ismissingfile = true ;
  end
end
