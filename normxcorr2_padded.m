function result = normxcorr2_padded(k, im, pad_style)
    % normalize k
    kc = k - mean(k, 'all') ;  % k centered
    Skk = sum(kc.^2, 'all') ;
    if Skk == 0 ,
        error('Kernel variance is zero')
    end
    kz = kc ./ sqrt(Skk) ;  % k, but z-scored
    [km, kn] = size(kz) ;
    ipad = ceil(km/2) ;
    jpad = ceil(kn/2) ;
    im_padded = padarray(im, [ipad jpad], pad_style) ;
    is_valid = logical(padarray(ones(size(im)), [ipad jpad])) ;  % true iff element in padded image should be kept when discarding pad
    result_padded_padded = normxcorr2(k, im_padded) ;  % this gives the 'full' correlution
    % Need to trim the 'full' result down to the 'same' result (same as im_padded,
    % of course)
    i_lo_full_pad = floor(km/2) ;
    i_hi_full_pad = floor((km-1)/2) ;
    j_lo_full_pad = floor(kn/2) ;
    j_hi_full_pad = floor((kn-1)/2) ;
    result_padded = result_padded_padded(i_lo_full_pad+1:end-i_hi_full_pad, j_lo_full_pad+1:end-j_hi_full_pad) ;  % same size as im_padded
    result = reshape(result_padded(is_valid), size(im)) ;    
end
