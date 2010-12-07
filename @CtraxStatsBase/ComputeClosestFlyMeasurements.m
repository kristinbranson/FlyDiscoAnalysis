% obj.ComputeClosestFlyMeasurements()
% Compute per-frame measurements of distances, angles to the nearest fly
% to each fly, as well as speeds related to these coordinate systems
% for all currently loaded trajectories. If new experiment directories
% are added after ComputeClosestFlyMeasurements is called, closest
% fly-derived measurements will be computed for the new trajectories
% when the experiment is added. The closest fly-based derived
% measurements will be stored to <expdir>/<closestflyfilestr>.

function ComputeClosestFlyMeasurements(obj)

if obj.didComputeClosestFlyMeasurements,
  return;
end

trx = [];
params = struct;
params.fov = obj.fov;
for n = 1:obj.nexpdirs,
  
  flies = obj.movie2flies{n};
  [trxnew,units] = ComputeClosestFlyMeasurements(obj.trx(flies),params,obj.closestflyfiles{n});
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

obj.didComputeClosestFlyMeasurements = true;
