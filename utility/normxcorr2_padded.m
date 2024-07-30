function result = normxcorr2_padded(k, im, pad_style)

% Make sure the kernel variance is nonzero.
% Then make the kernel zero-mean and unit-norm, why not.
kc = k - mean(k, 'all') ;  % k centered
Skk = sum(kc.^2, 'all') ;
if Skk == 0 ,
  error('Kernel variance is zero')
end
kz = kc ./ sqrt(Skk) ;  % k, but z-scored

% Pad im enough that the valid area is all of the original im.
[km, kn] = size(kz) ;
ipad = ceil(km/2) ;
jpad = ceil(kn/2) ;
im_padded = padarray(im, [ipad jpad], pad_style) ;

% Make a mark representing the non-pad parts of im_padded.
is_valid = logical(padarray(ones(size(im)), [ipad jpad])) ;

% Use the builtin normxcorr2() function to do the correlution on the padded
% image.  normxcorr2() always returns the 'full' correlution, so the result
% has ever more padding.
result_padded_padded = normxcorr2(k, im_padded) ;

% Need to trim off the padding added by normxcorr2().
% This is like trimming the results of conv2(x,y,'full') down to the results
% of conv2(x,y,'same').  Result will be same size as im_padded.
i_lo_full_pad = floor(km/2) ;
i_hi_full_pad = floor((km-1)/2) ;
j_lo_full_pad = floor(kn/2) ;
j_hi_full_pad = floor((kn-1)/2) ;
result_padded = result_padded_padded(i_lo_full_pad+1:end-i_hi_full_pad, j_lo_full_pad+1:end-j_hi_full_pad) ;  % same size as im_padded

% Finally, use the is_valid array to trim off the pad we added with
% padarray().
result = reshape(result_padded(is_valid), size(im)) ;

end  % function
