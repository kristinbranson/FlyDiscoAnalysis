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
  ns = obj.expdir2n(expdirs);
end

fprintf('Deleting the following speed measurment files:\n');
for n = ns,
  if ~exist(obj.speedfiles{n},'file'),
    continue;
  end
  fprintf('  %s\n',obj.speedfiles{n});
  delete(obj.speedfiles{n});
end

obj.didComputeClosestFlyMeasurements = false;