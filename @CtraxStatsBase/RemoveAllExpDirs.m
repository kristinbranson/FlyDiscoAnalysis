% obj.RemoveAllExpDirs()
% Removes all the data stored in obj, closes the relevant file handles.
function RemoveAllExpDirs(obj)

expdirs = obj.expdir_bases;

for n = 1:length(expdirs),
  obj.RemoveExpDir(expdirs{n});
end
  
