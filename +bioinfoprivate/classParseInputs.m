function obj = classParseInputs(obj, matlabfile, varargin)
%CLASSPARSEINPUTS Parse input PV pairs for classes. 
% 
%   CLASSPARSEINPUTS(OBJ, MATLABFILE) parses and validates input
%   property/value pair for object OBJ, and updates the object properties.

% Copyright 2009-2010 The MathWorks, Inc.


if nargin < 3
    return;
end

% Check for the right number of inputs
if rem(nargin-2, 2)== 1
    bioinfoprivate.bioclserror(matlabfile, matlabfile, 'IncorrectNumberOfArguments', ...
        ['Incorrect number of arguments to ' class(obj)])
end

obj = obj.set(varargin{:});
end % classParseInputs method
