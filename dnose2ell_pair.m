function d = dnose2ell_pair(trx,fly1,fly2)

% initialize
d = nan(1,trx(fly1).nframes);

% get start and end frames of overlap
t0 = max(trx(fly1).firstframe,trx(fly2).firstframe);
t1 = min(trx(fly1).endframe,trx(fly2).endframe);
  
% no overlap
if t1 < t0, 
  return;
end

% position of nose1
xnose = trx(fly1).x_mm + 2*trx(fly1).a_mm.*cos(trx(fly1).theta_mm);
ynose = trx(fly1).y_mm + 2*trx(fly1).a_mm.*sin(trx(fly1).theta_mm);

for t = t0:t1,
  i = t + trx(fly1).off;
  j = t + trx(fly2).off;
  d(i) = ellipsedist_hack(trx(fly2).x_mm(j),trx(fly2).y_mm(j),...
    trx(fly2).a_mm(j),trx(fly2).b_mm(j),trx(fly2).theta_mm(j),...
    xnose(i),ynose(i),50);
end