% obj.ComputeSpeedMeasurements()
% Compute per-frame measurements of speeds and accelerations for all
% currently loaded trajectories. If new experiment directories are
% added after ComputeSpeedMeasurements is called, speeds will
% be computed for the new trajectories when the experiment is added.
% The speed measurements will be stored to <expdir>/<speedfilestr>.

function ComputeSpeedMeasurements(obj)

if obj.didComputeSpeedMeasurements,
  return;
end

trx = [];
params = struct;
params.thetafil = obj.thetafil;
for n = 1:obj.nexpdirs,
  
  flies = obj.movie2flies{n};
  [trxnew,units] = ComputeSpeedMeasurements(obj.trx(flies),params,obj.speedfiles{n});
  trx = structappend(trx,trxnew);
  
end

obj.trx = trx;
clear trx;

% add new units
fns = fieldnames(units);
for i = 1:length(fns),
  fn = fns{i};
  obj.units.(fn) = units.(fn);
end

obj.didComputeSpeedMeasurements = true;
