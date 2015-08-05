function pvPairs = validateClassPVNames(obj, errmfilename, varargin)
%VALIDATECLASSPVPAIR Validate parameter names match class property names.
%
%   PVPAIR = VALIDATECLASSPVPAIR(OBJ, ERRFILENAME, INPUTPVPAIR) Validates
%   the input parameter names of parameter/value pair matches OBJ class
%   property names. The input parameter names can be case insensitive. The
%   returned PVPAIR parameter names are replaced with corresponding class
%   property names. The PVPAIR can be pass into the set method of the
%   class.  
%
%   Note: This function is useful for classes using hgsetget. 

% Copyright 2009 The MathWorks, Inc.


% Check for the right number of inputs
if rem(nargin-2,2)== 1
    bioinfoprivate.bioerror(errmfilename,...
        'IncorrectNumberOfArguments', ...
        'Incorrect number of arguments to %s.', errmfilename)
end

% Allowed inputs
propertyNames = fieldnames(obj);
errMsg = '';
pvPairs = varargin;
for j=1:2:nargin-2
    pName = varargin{j};
    
    if ~(ischar(pName) && isvector(pName) && (size(pName,1)==1))
        bioinfoprivate.bioclserror(errmfilename, errmfilename,...
            'InvalidPropertyName', 'Invalid property name.');
    end
    
    % Using lower to be case insensitive
    k = find(nominal(lower(propertyNames)) == lower(pName)); 
    if isempty(k) || length(k)>1
        errMsg = [errMsg pName ',']; %#ok
    else
        pvPairs{j} = propertyNames{k};
    end
end
if ~isempty(errMsg)
    bioinfoprivate.bioclserror(errmfilename, errmfilename,...
        'InvalidParameterName', ...
        'Unknown or ambiguous parameter name(s) %s.', errMsg(1:end-1))
end
end
