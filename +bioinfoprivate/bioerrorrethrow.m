function bioerrorrethrow(matlabfile, orignME, varargin)
%BIOERRORRETHROW rethrows error as caller for functions. 
% 
%   BIOERRORRETHROW(MATLABFILE, ORIGNME) rethrows MException with error ID
%   and error MSG for MATLABFILE from original MException ORIGME, but the
%   prefix of error identifier is replaced by: bioinfo:filename.

% Copyright 2008-2010 The MathWorks, Inc.


idx = strfind(orignME.identifier, ':');
idx = idx(end) + 1;
    
x = bioinfoprivate.bioexception(matlabfile, orignME.identifier(idx:end),...
        orignME.message, varargin{:});
    
x.throwAsCaller;
end % error method
