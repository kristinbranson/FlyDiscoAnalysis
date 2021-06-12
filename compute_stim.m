function [data,units] = compute_stim(trx,n,stimtype,stimnums)

ind = false(1,trx.movie_nframes(n));
for i = 1:numel(stimnums),
  if strcmpi(stimtype,'on'),
    t0 = trx.indicatorLED{n}.starton(i);
    t1 = trx.indicatorLED{n}.endon(i);
  else
    t0 = trx.indicatorLED{n}.startoff(i);
    t1 = trx.indicatorLED{n}.endoff(i);
  end
  ind(t0:t1) = true;
end
data = cell(1,trx.nfliespermovie(n));
flies = trx.exp2flies{n};
for flyi = 1:numel(flies),
  fly = flies(flyi);
  data{flyi} = ind(trx.firstframes(fly):trx.endframes(fly));
end
units = parseunits('unit');