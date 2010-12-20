% [ns,flies] = obj.IntersectFliesExpdirs(flies,expdirs)
% Returns the flies, experiments that are both within flies and
% experiments
function ns = expdir2n(obj,expdirs)

[didfind,ns] = ismember(expdirs,obj.expdir_bases);
if any(~didfind),
  warning(['The following experiments are not loaded:\n',sprintf('%s\n',expdirs{~didfind})]);
end