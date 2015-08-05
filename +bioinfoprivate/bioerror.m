function bioerror(matlabfile, id, msg,varargin)
%BIOERROR throws error as caller
% 
%   BIOERROR(MATLABFILE, ID, MSG) throws MException with error ID and
%   error MSG for MATLABFILE.

% Copyright 2008-2010 The MathWorks, Inc.


x = bioinfoprivate.bioexception(matlabfile, id, msg, varargin{:});
x.throwAsCaller;
end % error method
