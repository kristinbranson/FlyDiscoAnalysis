function [k, pval] = optPartialMatch(pval, oktypes, parametername, matlabfile)
%OPTPARTIALMATCH Helper function that looks for partial matches of string
% options in a list of inputs and returns the complete string for the
% selected option.
%
%   [K, PVAL] = OPTPARTIALMATCH(PVAL, OKTYPES, PNAME, MATLABFILE) given
%   an input string PVAL and a cell with possible matches OKTYPES, finds
%   the matching unambiguous string and returns K, the index of the match
%   and the complete string for the selected option. Function returns an
%   error in the case of an ambiguous match or an failed match. PNAME and
%   MATLABFILE are used to build the error ID.

% Copyright 2010 The MathWorks, Inc.


k = find(strncmpi(pval, oktypes, numel(pval)));
if numel(k) == 1
    pval = oktypes{k};
    return
end

% Prepare the error message
if numel(oktypes)==1
   vostr = sprintf('Valid option is ''%s''.',oktypes{1});
elseif numel(oktypes)==2
   vostr = sprintf('Valid options are ''%s'' or ''%s''.',oktypes{1},oktypes{2});
else
   vostr = ['Valid options are ', ...
            sprintf('''%s'', ',oktypes{1:end-1}), ...
            sprintf('or ''%s''.',oktypes{end})];
end

% Throw the respective exception
if isempty(k)
    x = bioinfoprivate.bioexception(matlabfile, ['Unknown',parametername],...
        'Unknown option for %s. %s',upper(parametername),vostr);
    x.throwAsCaller;
elseif length(k)>1
    x = bioinfoprivate.bioexception(matlabfile, ['Ambiguous',parametername],...
        'Ambiguous option for %s. %s',upper(parametername),vostr);
    x.throwAsCaller;
end

end % optPartialMatch method
