function x = bioexception(matlabfile, id, msg, varargin)
%BIOEXCEPTION returns MException
%
%   X = BIOEXCEPTION(MATLABFILE, ID, MSG) returns MException X with identifier
%   ID and message MSG for MATLABFILE.

% Copyright 2008-2010 The MathWorks, Inc.


msgId = sprintf('bioinfo:%s:%s', matlabfile, id);
x = MException(msgId, msg, varargin{:});
end % bioexception
