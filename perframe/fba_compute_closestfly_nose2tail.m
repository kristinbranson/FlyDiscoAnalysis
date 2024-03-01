% closest fly, based on dnose2tail
function [data,units,mind] = fba_compute_closestfly_nose2tail(trx,n,dosave_d)

if nargin < 3,
  dosave_d = true;
end

flies = trx.exp2flies{n};
nflies = numel(flies);
closestfly = cell(1,nflies);
mind = cell(1,nflies);

for i1 = 1:nflies,
  fly1 = flies(i1);
  fprintf('fly1 = %d\n',fly1);
  fly1_chamber_index = trx(fly1).chamber_index ;  % scalar
  d = nan(nflies,trx(fly1).nframes);
  for i2 = 1:nflies,
    if i1 == i2,
      continue;
    end
    fly2 = flies(i2);
    fly2_chamber_index = trx(fly2).chamber_index ;  % scalar
    if fly2_chamber_index ~= fly1_chamber_index ,
      continue
    end
    d(i2,:) = dnose2tail_pair(trx,fly1,fly2);
  end
  [mind{i1},closesti] = min(d,[],1);
  closestfly{i1} = flies(closesti);
  closestfly{i1}(isnan(mind{i1})) = nan;
end

% so that we don't compute dcenter twice
if dosave_d,
  data = mind;
  units = parseunits('mm');
  filename = trx.GetPerFrameFile('dnose2tail',n);
  save(filename,'data','units');
end

data = closestfly;
units = parseunits('unit');
