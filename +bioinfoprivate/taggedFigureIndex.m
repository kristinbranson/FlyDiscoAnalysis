function n = taggedFigureIndex(tag, baseName)
%TAGGEDFIGUREINDEX Returns the lowest number for a figure with a tag.
% 
%   I = TAGGEDFIGUREINDEX(TAG, BASENAME) computes the lowest integer number
%   I for a figure with a specified tag TAG and a name BASENAME. The index
%   I can be added to the BASENAME for the next figure to keep count of the
%   figures with the same tag. For example, Clustergram 1 and Clustergram 2
%   are the figure already opened with the tag of clustergram, to open
%   another clustergram, the new figure name would be Clustergram 3. 

% Copyright 2009 The MathWorks, Inc.


%== Find all the used numbers so far
allFigs = findall(0,'tag', tag);
usedNumbers = zeros(1,numel(allFigs)+1);
baseName = [strtrim(baseName) ' '];
baseLen = length(baseName);
for i = 1:numel(allFigs)
    str = get(allFigs(i),'Name');
    usedNumbers(i) = str2double(str(baseLen:end));
end
%== Find the next index
% The rule is that we find the lowest integer value (non-zero and positive)
% not yet prescribed to a figure tagged with TAG. This is the same way
% MATLAB figures behave. 
n = min(setdiff(1:(max(usedNumbers)+1),usedNumbers));
end
