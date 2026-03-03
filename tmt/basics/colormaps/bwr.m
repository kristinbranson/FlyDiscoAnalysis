function cmap=bwr(n)

if ~exist('n', 'var') || isempty(n) ,
    n = 256 ;
end

% n is the number of colors
clr_pos=[1 0 0];
clr_zero = [1 1 1] ;
clr_neg=[0 0 1];

clr_pos_lab=srgb2lab(clr_pos);
clr_zero_lab=srgb2lab(clr_zero);
clr_neg_lab=srgb2lab(clr_neg);

n_half=ceil(n/2);
scale=linspace(0,1,n_half)';
cmap_pos_lab = scale .* repmat(clr_pos_lab,[n_half 1]) + (1-scale) .* repmat(clr_zero_lab,[n_half 1]) ;
cmap_neg_lab = scale .* repmat(clr_neg_lab,[n_half 1]) + (1-scale) .* repmat(clr_zero_lab,[n_half 1]) ;
if mod(n,2)==0 ,
    cmap_lab=[flipud(cmap_neg_lab);cmap_pos_lab];
else
    cmap_lab=[flipud(cmap_neg_lab); clr_zero_lab; cmap_pos_lab];
end
cmap=lab2srgb(cmap_lab);
cmap=max(min(cmap,1),0);
