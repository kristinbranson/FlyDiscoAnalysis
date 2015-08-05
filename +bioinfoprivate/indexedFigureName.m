function fname = indexedFigureName(tag, basename)
%INDEXEDFIGURENAME Returns the name with lowest number for a figure.
% 
%   FIGNAME = INDEXEDFIGURENAME(TAG, BASENAME) returns a name FIGNAME with
%   the lowest integer number a figure with a specified tag TAG and a name
%   BASENAME. The index is added after the BASENAME for the next figure to
%   keep count of the figures with the same tag. For example, Clustergram
%   1, Clustergram 2 are the figure already opened with the tag of
%   clustergram, to open another clustergram, the new figure name would be
%   Clustergram 3. 

% Copyright 2009 The MathWorks, Inc.


%== Find all the used numbers so far
fname = sprintf([basename ' %d'], bioinfoprivate.taggedFigureIndex(tag, basename));
end 
