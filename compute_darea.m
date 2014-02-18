% change in area
function [data,units] = compute_darea(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  data{i} = diff(trx(fly).area,1,2)./trx(fly).dt;
end
units = parseunits('mm^2/s');

