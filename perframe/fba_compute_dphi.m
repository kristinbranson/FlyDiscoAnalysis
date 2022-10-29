% change in velocity direction
function [data,units] = compute_dphi(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  diffphi = diff(trx(fly).phi,1,2);
  flydt = trx(fly).dt;
  data{i} = modrange(diffphi,-pi,pi)./flydt;
end
units = parseunits('rad/s');

