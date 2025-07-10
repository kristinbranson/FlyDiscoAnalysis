function RotateFliesInTime(trx)

nflies = numel(trx);
firstframe = min([trx.firstframe]);
endframe = max([trx.endframe]);
nframes = endframe-firstframe+1;

% how many flies are there at any one time
nfliesperframe = zeros(1,nframes);
for fly = 1:nflies,
  nfliesperframe(trx(fly).firstframe:trx(fly).endframe) = ...
    nfliesperframe(trx(fly).firstframe:trx(fly).endframe) + 1;
end