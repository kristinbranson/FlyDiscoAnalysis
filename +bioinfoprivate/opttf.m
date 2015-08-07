function tf = opttf(pval, okarg, matlabfile)
%OPTTF determines whether input options are true or false
%
%   TF = OPTTF(PVAL, OKARG, MATLABFILE) evaluates the value PVAL of
%   property OKARG for logical TF and errors if the intended valus is
%   invalid. 

% Copyright 2008-2010 The MathWorks, Inc.


if islogical(pval)
    tf = all(pval);
    return
end
if isnumeric(pval)
    tf = all(pval~=0);
    return
end
if ischar(pval)
    truevals = {'true','yes','on','t'};
    k = any(strcmpi(pval,truevals));
    if k
        tf = true;
        return
    end
    falsevals = {'false','no','off','f'};
    k = any(strcmpi(pval,falsevals));
    if k
        tf = false;
        return
    end
end
if nargin == 1
    % return empty if unknown value
    tf = logical([]);
else
    okarg(1) = upper(okarg(1));
    x = bioinfoprivate.bioexception(matlabfile, [okarg 'OptionNotLogical'],...
        '%s must be a logical value, true or false.', upper(okarg));
    x.throwAsCaller;
end
end % opttf function
