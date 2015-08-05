function [name, k] = validateClassPropName(obj, name, matlabfile, varargin)
%VALIDATECLASSPROPNAME Validate class property name. 
% 
%   [NAME, K] = VALIDATECLASSPROPNAME(OBJ, NAME, MATLABFILE) validates NAME
%   matches the property names of an object OBJ, and returns the validated
%   name and index.
%
%   [NAME, K] = VALIDATECLASSPROPNAME(OBJ, NAME, MATLABFILE, EXTRANAMES)
%   validates EXTRANAMES.

% Copyright 2009-2010 The MathWorks, Inc.


% This matches names against this list of property names.
clsname = class(obj);
dotIdx = strfind(clsname, '.');
clsname = clsname(dotIdx(end)+1:end);
if ~(ischar(name) && isvector(name) && (size(name,1)==1))
    bioinfoprivate.bioclserror(clsname, matlabfile, 'InvalidPropertyName', ...
          'Invalid property name.');
end

% Add extra names in varargin
if isempty(varargin)
    propertyNames = fieldnames(obj);
else
    extra = varargin{:};
    propertyNames = [fieldnames(obj); extra(:)];
end
    
k = bioinfoprivate.pvpair(name, [], propertyNames, [clsname ':' matlabfile]);
name = propertyNames{k};
end 
