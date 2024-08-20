function result = tweak_rgb_image_size(im, desired_height, desired_width)
    % result is like im, but of size [desired_height desired_width 3], 
    % achieved by padding/trimming im in x and y
    
    original_height = size(im,1) ;
    original_width = size(im,2) ;
    class_name = class(im) ;
    
    % trim/pad in y
    if original_height < desired_height ,
        dheight_lo = floor((desired_height-original_height)/2) ;
        dheight_hi = (desired_height-original_height)-dheight_lo ;
        height_lo_pad = zeros([dheight_lo,original_width,3],class_name) ;
        height_hi_pad = zeros([dheight_hi,original_width,3],class_name) ;
        pre_result = cat(1, height_lo_pad, im, height_hi_pad) ;
    elseif original_height > desired_height ,
        dheight_lo = floor((original_height-desired_height)/2) ;
        dheight_hi = (original_height-desired_height)-dheight_lo ;
        pre_result = im(1+dheight_lo:end-dheight_hi, :, :) ;
    else
        pre_result = im ;
    end
    
    % trim/pad in x
    if original_width < desired_width ,
        dwidth_lo = floor((desired_width-original_width)/2) ;
        dwidth_hi = (desired_width-original_width)-dwidth_lo ;
        width_lo_pad = zeros([desired_height,dwidth_lo,3],class_name) ;
        width_hi_pad = zeros([desired_height,dwidth_hi,3],class_name) ;
        result = cat(2, width_lo_pad, pre_result, width_hi_pad) ;
    elseif original_width > desired_width ,
        dwidth_lo = floor((original_width-desired_width)/2) ;
        dwidth_hi = (original_width-desired_width)-dwidth_lo ;
        result = pre_result(:, 1+dwidth_lo:end-dwidth_hi, :) ;
    else
        result = pre_result ;
    end
end
