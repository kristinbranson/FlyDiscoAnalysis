% obj.DeleteLandmarkMeasurementFiles([expdirs])
% Deletes the files created during the ComputeLandmarkMeasurements step
% -- the file containing the landmark-based derived measurements
% <expdir>/<landmarksfilestr>. If expdirs is input, then the
% landmark files for the specified expdirs will be deleted. If not
% specified, then landmark files for all loaded experiments will be
% deleted.

function DeleteLandmarkMeasurementFiles(obj,expdirs)

if ~exist('expdirs','var'),
  ns = 1:obj.nexpdirs;
else
  [didfind,ns] = ismember(expdirs,obj.expdirs);
  for i = find(~didfind),
    warning('Expdir %s not loaded.\n',expdirs{i});
  end
  ns = ns(didfind);
end

fprintf('Deleting the following landmark measurment files:\n');
for n = ns,
  if ~exist(obj.landmarksfiles{n},'file'),
    continue;
  end
  fprintf('  %s\n',obj.landmarksfiles{n});
  delete(obj.landmarksfiles{n});
end