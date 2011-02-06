% magnitude of difference in velocity between fly and closest fly based on
% type. 
function [data,units] = compute_magveldiff(trx,n,type)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);

for i1 = 1:nflies,
  fly1 = flies(i1);
  closestfly = trx(fly1).(['closestfly_',type]);
  for i2 = 1:nflies,
    fly2 = flies(i2);
    if i1 == i2, continue; end
    idx = find(closestfly == fly2);
    if isempty(idx), continue; end
    off = trx(fly1).firstframe - trx(fly2).firstframe;
    dx = trx(fly2).x_mm(off+idx) - trx(fly1).x_mm(idx);
    dy = trx(fly2).y_mm(off+idx) - trx(fly1).y_mm(idx);
    data{i1}(idx) = sqrt(dx.^2 + dy.^2);
  end
end
data{i1} = data{i1} ./ trx(fly1).dt;
units = parseunits('mm/s');