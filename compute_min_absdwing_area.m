function [data,units] = compute_min_absdwing_area(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  data{i} = min(abs(diff(trx(fly).wing_areal_mm,1,2)),abs(diff(trx(fly).wing_arear_mm,1,2))) ./ trx(fly).dt;
  
end
units = parseunits('mm^2/s');

