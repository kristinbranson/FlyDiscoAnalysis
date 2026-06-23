% smoothed magnitude of velocity of center
% Speed signal for speed-conditional binning of locomotion metrics: a
% Gaussian-weighted moving average of velmag_ctr (Wilson 2024 idiom).
% Window hardcoded at 15 frames (=100ms @150fps), the value validated in
% Phase 0 against the control velmag distribution. Smooths the velmag_ctr
% signal directly; co-indexed with velmag_ctr (length nframes-1).
function [data,units] = compute_smoothvelmag_ctr(trx,n)

smoothwin = 15;  % frames; validated speed-bin smoothing window (Phase 0)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  vc = trx(fly).velmag_ctr;
  if isempty(vc),
    data{i} = [];
  else
    data{i} = smoothdata(vc,2,'gaussian',smoothwin,'omitnan');
  end
end
units = parseunits('mm/s');
