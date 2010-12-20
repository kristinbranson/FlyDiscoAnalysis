% ComputeLandmarkMeasurements(trx,params,landmarks_file)
% Compute per-frame measurements of distances, angles to landmarks in
% the arena
function [trx,units] = ComputeLandmarkMeasurements(trx,params,landmarks_file)

nflies = length(trx);
fns = {'dist2wall','wallangle','angle2wall','arena_r','ddist2wall','dwallangle','dangle2wall'};

isfile = exist('landmarks_file','var');

if isfile && exist(landmarks_file,'file'),
  fprintf('Loading landmark measurements from file %s\n',landmarks_file);
  load(landmarks_file,'landmarks','units');
else

  fprintf('Computing landmark measurements from file %s\n',landmarks_file);
  landmarks = structallocate(fns,[1,nflies]);
  units = struct;
  units.dist2wall = parseunits('mm');
  units.wallangle = parseunits('rad');
  units.angle2wall = parseunits('rad');

  for fly = 1:nflies,
    
    % distance to wall
    dcenter = sqrt((trx(fly).x_mm - params.arena_center_mm(1)).^2 + ...
      (trx(fly).y_mm - params.arena_center_mm(2)).^2);
    landmarks(fly).arena_r = dcenter;
    landmarks(fly).dist2wall = params.arena_radius_mm - dcenter;
    
    % polar angle of closest point on the wall
    landmarks(fly).wallangle = atan2(trx(fly).y_mm-params.arena_center_mm(2),...
      trx(fly).x_mm-params.arena_center_mm(1));
    
    % angle to closest point on the wall in the fly's coordinate system
    landmarks(fly).angle2wall = modrange(landmarks(fly).wallangle-trx(fly).theta_mm,-pi,pi);
    
    % change in distance to wall
    landmarks(fly).ddist2wall = diff(landmarks(fly).dist2wall)./trx(fly).dt;
    landmarks(fly).units.ddist2wall = parseunits('mm/s');

    % change in polar angle to closest point on wall
    landmarks(fly).absdwallangle = ...
      abs(modrange(diff(landmarks(fly).wallangle),-pi,pi))./trx(fly).dt;
    landmarks(fly).units.absdwallangle = parseunits('rad/s');

    % change in angle to closest point on wall in fly's coordinate system
    landmarks(fly).dangle2wall = ...
      modrange(diff(landmarks(fly).angle2wall),-pi,pi)./trx(fly).dt;
    landmarks(fly).units.absdangle2wall = parseunits('rad/s');
    
  end

  if isfile,
    fprintf('Saving landmark measurements from file %s\n',landmarks_file);
    save(landmarks_file,'landmarks','units');
  end
end

% add these fields to trx
for i = 1:length(fns),
  fn = fns{i};
  for fly = 1:nflies,
    trx(fly).(fn) = landmarks(fly).(fn);
  end
end
