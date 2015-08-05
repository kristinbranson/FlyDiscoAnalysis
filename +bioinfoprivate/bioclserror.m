function bioclserror(clsname, matlabfile, id, msg, varargin)
%BIOCLSERROR throws error as caller for functions in a MCOS class directory. 
% 
%   BIOCLSERROR(CLSNAME, MATLABFILE, ID, MSG) throws MException with error
%   ID and error MSG for MATLABFILE in class CLSNAME.

% Copyright 2008-2010 The MathWorks, Inc.


x = bioinfoprivate.bioexception([clsname ':' matlabfile], id, msg, varargin{:});
x.throwAsCaller;
end % error method
