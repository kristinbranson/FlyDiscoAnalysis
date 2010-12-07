% obj.ComputeLandmarkMeasurements()
% Compute per-frame measurements of distances, angles to landmarks in
% the arena, as well as speeds related to these coordinate systems for
% all currently loaded trajectories. If new experiment directories are
% added after ComputeLandmarkMeasurements is called, landmark-derived
% measurements will be computed for the new trajectories when the
% experiment is added. Currently, the only landmark is the circular
% arena wall, which is assumed to have an origin arena_center_mm and
% radius arena_radius_mm in the registered coordinate system. The
% landmark-based derived measurements will be stored to
% <expdir>/<landmarksfilestr>.

function ComputeLandmarkMeasurements(obj)

if obj.didComputeLandmarkMeasurements,
  return;
end

trx = [];
params = struct;
params.arena_center_mm = obj.arena_center_mm;
params.arena_radius_mm = obj.arena_radius_mm;
for n = 1:obj.nexpdirs,
  
  flies = obj.movie2flies{n};
  [trxnew,units] = ComputeLandmarkMeasurements(obj.trx(flies),params,obj.landmarksfiles{n});
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

obj.didComputeLandmarkMeasurements = true;
