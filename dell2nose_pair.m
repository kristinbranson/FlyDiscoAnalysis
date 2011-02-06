function d = dell2nose_pair(trx,fly1,fly2)

% initialize
d = nan(1,trx(fly1).nframes);

% get start and end frames of overlap
t0 = max(trx(fly1).firstframe,trx(fly2).firstframe);
t1 = min(trx(fly1).endframe,trx(fly2).endframe);
  
% no overlap
if t1 < t0, 
  return;
end

% position of nose2
xnose = trx(fly2).x_mm + 2*trx(fly1).a_mm.*cos(trx(fly2).theta_mm);
ynose = trx(fly2).y_mm + 2*trx(fly1).a_mm.*sin(trx(fly2).theta_mm);

for t = t0:t1,
  i = t + trx(fly1).off;
  j = t + trx(fly2).off;
  d(i) = ellipsedist_hack(trx(fly1).x_mm(i),trx(fly1).y_mm(i),...
    trx(fly1).a_mm(i),trx(fly1).b_mm(i),trx(fly1).theta_mm(i),...
    xnose(j),ynose(j),50);
end