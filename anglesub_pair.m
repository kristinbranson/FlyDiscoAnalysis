function anglesub = anglesub_pair(trx,fly1,fly2)

% initialize
anglesub = nan(1,trx(fly1).nframes);

% get start and end frames of overlap
t0 = max(trx(fly1).firstframe,trx(fly2).firstframe);
t1 = min(trx(fly1).endframe,trx(fly2).endframe);
  
% no overlap
if t1 < t0, 
  return;
end

for t = t0:t1,
  i = t + trx(fly1).off;
  j = t + trx(fly2).off;
  anglesub(i) = anglesubtended(...
    trx(fly1).x_mm1(i),trx(fly1).y_mm1(i),trx(fly1).a_mm1(i),trx(fly1).b_mm1(i),trx(fly1).theta_mm1(i),...
    trx(fly2).x_mm2(j),trx(fly2).y_mm2(j),trx(fly2).a_mm2(j),trx(fly2).b_mm2(j),trx(fly2).theta_mm2(j),...
    trx.perframe_params.fov);
end