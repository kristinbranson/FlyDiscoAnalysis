function [i, j, m] = imargmax(im)
    % Get the (i,j) coords at which the max of image im occurs.  Also returns the
    % max as third return value.
    [m, ij_serial] = max(im, [], 'all') ;
    [i, j] = ind2sub(size(im), ij_serial) ;
end
