im = zeros(500,500) ;
im = max(im + normrnd(0,0.01,size(im)), 0) ;
i_target = 100 
j_target = 50 
im(i_target,j_target) = 1 ;
k = zeros(11,11) ;
k(6,6) = 1 ;
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
