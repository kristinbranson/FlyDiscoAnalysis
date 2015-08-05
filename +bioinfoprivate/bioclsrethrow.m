function bioclsrethrow(clsname, matlabfile, orignME, varargin)
%BIOCLSRETHROW rethrows error as caller for functions in a MCOS class
%   directory. 
% 
%   BIOCLSRETHROW(CLSNAME, MATLABFILE, ORIGNME) rethrows MException with
%   error ID and error MSG for MATLABFILE in class CLSNAME from OLDEXCEPTION,
%   but the prefix of error identifier is replaced by:
%   bioinfo:classname:filename.

% Copyright 2008-2010 The MathWorks, Inc.


idx = strfind(orignME.identifier, ':');
idx = idx(end) + 1;
    
x = bioinfoprivate.bioexception([clsname ':' matlabfile], orignME.identifier(idx:end),...
        orignME.message, varargin{:});
    
x.throwAsCaller;
end % error method
