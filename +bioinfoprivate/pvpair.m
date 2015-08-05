function [k, pval] = pvpair(pname, theVal, okargs, mfile)
%PVPAIR Helper function that looks for partial matches of parameter
% names in a list of inputs and returns the parameter/value pair and
% matching number.
%
%   [K, PVAL] = PVPAIR(PNAME, THEVAL, OKARGS) given input string PNAME,
%   and corresponding value, THEVAL, finds matching name in the OKARGS
%   list. Returns K, the index of the match, and PVAL, the parameter
%   value.

% Copyright 2008 The MathWorks, Inc.


if ~ischar(pname) || ~isrow(pname)
    x = bioinfoprivate.bioexception(mfile, 'InvalidParameterName',...
        'Parameter name must be a character string.');
    x.throwAsCaller;
end

k = find(strncmpi(pname, okargs,numel(pname)));

if numel(k) == 1
    pval = theVal;
    return
end

if isempty(k)
    x = bioinfoprivate.bioexception(mfile, 'UnknownParameterName',...
        'Unknown parameter name: %s.',pname);
    x.throwAsCaller;

elseif length(k)>1
    x = bioinfoprivate.bioexception(mfile, 'AmbiguousParameterName',...
        'Ambiguous parameter name: %s.',pname);
    x.throwAsCaller;
end
end %pvpair method
