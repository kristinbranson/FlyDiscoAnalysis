function [isfore,diffim,l] = BackSub(im,ann)

BG_TYPE_LIGHTONDARK = 0;
BG_TYPE_DARKONLIGHT = 1;
BG_TYPE_OTHER = 2;

switch ann.bg_type
  case BG_TYPE_LIGHTONDARK,
    diffim = max(0,im - ann.background_center);
  case BG_TYPE_DARKONLIGHT,
    diffim = max(0,ann.background_center - im);
  case BG_TYPE_OTHER,
    diffim = abs(ann.background_center - im);
end

diffim = diffim ./ ann.background_dev;
isfore = diffim >= ann.n_bg_std_thresh_low;
isfore_high = diffim >= ann.n_bg_std_thresh;
[l,n] = bwlabel(isfore);
for i = 1:n,
  if ~any(isfore_high(l == i)),
    isfore(l==i) = false;
    l(l==i) = 0;
  end
end
[~,~,l] = unique(l);
l = reshape(l,size(diffim)) - 1;