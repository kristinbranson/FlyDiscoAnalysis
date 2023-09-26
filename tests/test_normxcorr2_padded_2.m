im = exprnd(100,[500 500]) ;
i_target = 100 
j_target = 50 
radius = 5 
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
