% obj.DeleteSpeedMeasurementFiles([expdirs])
% Deletes the files created during the ComputeSpeedMeasurements step --
% the file containing the speed measurements <expdir>/<speedfilestr>.
% If expdirs is input, then the speed files for the specified expdirs
% will be deleted. If not specified, then speed files for all loaded
% experiments will be deleted.

function DeleteSpeedMeasurementFiles(obj,expdirs)

if ~exist('expdirs','var'),
  ns = 1:obj.nexpdirs;
else
  [didfind,ns] = ismember(expdirs,obj.expdirs);
  for i = find(~didfind),
    warning('Expdir %s not loaded.\n',expdirs{i});
  end
  ns = ns(didfind);
end

fprintf('Deleting the following speed measurment files:\n');
for n = ns,
  if ~exist(obj.speedfiles{n},'file'),
    continue;
  end
  fprintf('  %s\n',obj.speedfiles{n});
  delete(obj.speedfiles{n});
end