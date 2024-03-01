function [data,units] = fba_compute_nflies_close(trx,n,nbodylengths_near)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);

if nargin < 3,
  nbodylengths_near = trx.perframe_params.nbodylengths_near;
end


for i1 = 1:nflies,
  fly1 = flies(i1);
  %fprintf('fly1 = %d\n',fly1);
  fly1_chamber_index = trx(fly1).chamber_index ;  % scalar
  data{i1} = zeros(1,trx(fly1).nframes);
  for i2 = 1:nflies,
    if i1 == i2,
      continue
    end
    fly2 = flies(i2);
    fly2_chamber_index = trx(fly2).chamber_index ;  % scalar
    if fly2_chamber_index ~= fly1_chamber_index ,
      continue
    end    
    idx = isclose_pair(trx,fly1,fly2,nbodylengths_near);
    data{i1}(idx) = data{i1}(idx) + 1;
  end
end

units = parseunits('unit');
