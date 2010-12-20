% [ns,flies] = obj.IntersectFliesExpdirs(flies,expdirs)
% Returns the flies, experiments that are both within flies and
% experiments
function [ns,flies] = IntersectFliesExpdirs(obj,flies,expdirs)

% take the intersection of specified flies and expdirs

% which expdirs will we look at?
ns = obj.expdir2n(expdirs);

% which flies does this correspond to?
allflies_perexp = [obj.movie2flies{ns}];
% take the intersection of the specified flies and the flies in the
% specified experiments
flies = intersect(flies,allflies_perexp);
% which experiments do these flies correspond to
ns = unique(obj.fly2movie(flies));
