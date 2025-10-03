% absolute value of sideways motion of center of rotation
function [data,units] = compute_absdu_ctr(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  
  data{i} = abs(trx(fly).du_ctr);

end
units = parseunits('mm/s');

