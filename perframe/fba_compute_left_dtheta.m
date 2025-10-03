% left sideways motion of the body
function [data,units] = compute_left_vel(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  
  data{i} = max(0,trx(fly).dv_ctr);

end
units = parseunits('mm/s');

