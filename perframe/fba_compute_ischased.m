function [data,units] = compute_ischased(trx,expi)

units = parseunits('unit');
flies = trx.exp2flies{expi};
nflies = numel(flies);
data = cell(1,nflies);

for fly1i = 1:nflies,
  fly1 = flies(fly1i);
  ischased = false(1,trx(fly1).nframes);
  for fly2 = 1:nflies,
    if fly1 == fly2,
      continue;
    end
    t0 = max(trx(fly1).firstframe,trx(fly2).firstframe);
    t1 = min(trx(fly1).endframe,trx(fly2).endframe);
    if t0 > t1,
      continue;
    end
    i0 = t0 - trx(fly1).firstframe + 1;
    i1 = t1 -trx(fly1).firstframe + 1;
    j0 = t0 - trx(fly2).firstframe + 1;
    j1 = t1 -trx(fly2).firstframe + 1;
    
    ischased_curr = trx(fly2).closestfly_nose2ell_angle_min30to30(j0:j1) == fly1i & ...
      trx(fly2).chase_labels(j0:j1);
    ischased(i0:i1) = ischased(i0:i1)|ischased_curr;
    
  end
  data{fly1} = ischased;
end

