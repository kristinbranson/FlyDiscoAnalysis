% obj.DeleteClosestFlyMeasurementFiles([expdirs])
% Deletes the files created during the ComputeClosestFlyMeasurements step
% -- the file containing the closest fly-based derived measurements
% <expdir>/<closestflyfilestr>. If expdirs is input, then the
% closest-fly files for the specified expdirs will be deleted. If not
% specified, then closest fly files for all loaded experiments will be
% deleted.

function DeleteClosestFlyMeasurementFiles(obj,expdirs)

if ~exist('expdirs','var'),
  ns = 1:obj.nexpdirs;
else
  [didfind,ns] = ismember(expdirs,obj.expdirs);
  for i = find(~didfind),
    warning('Expdir %s not loaded.\n',expdirs{i});
  end
  ns = ns(didfind);
end

fprintf('Deleting the following closestfly measurment files:\n');
for n = ns,
  if ~exist(obj.closestflyfiles{n},'file'),
    continue;
  end
  fprintf('  %s\n',obj.closestflyfiles{n});
  delete(obj.closestflyfiles{n});
end