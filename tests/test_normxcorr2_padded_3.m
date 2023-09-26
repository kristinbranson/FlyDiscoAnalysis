m = 500 
n = 500
i_target = 100 
j_target = 50
ayes = repmat((1:m)', [1 n]) ;
jays = repmat((1:n) , [m 1]) ;
sigma = 3 ;
radius = sigma*5 ;
im_pure =exp (-0.5 * ((ayes-i_target).^2 + (jays-j_target).^2) / sigma^2 ) ;
im = exprnd(0.1,[m n]) + im_pure ;
k = im(i_target-radius:i_target+radius, j_target-radius:j_target+radius) ;
result = normxcorr2_padded(k, im, 'replicate') ;

imglance(im)
title('im') ;
imglance(k)
title('k') ;
imglance(result)
title('result') ;
[result_max, result_argmax_serial] = max(result, [], 'all') ;
result_max
[i_result_argmax, j_result_argmax] = ind2sub(size(result), result_argmax_serial)
