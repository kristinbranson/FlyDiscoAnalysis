function [f,a,ih] = imglance(im, scale)
    % Get a quick look at an image, usually while debugging.
    % Behaves somewhat like imshow(), but imshow() is suboptimal in several ways,
    % IMHO.
    % Note to self: The purpose of this is to plot something reasonable in all cases.
    % Don't change it just becuse it doesn't plot something perfect in some weird case.
    if ~exist('scale', 'var') || isempty(scale) ,
        scale = '' ;  % Let imglance() decide about scaling
    end
    if strcmp(scale, 'scale') || strcmp(scale, 'dont') || strcmp(scale, '') ,
        % do nothing, all is well
    else
        error('scale_string should be one of: ''scale'', ''dont'' or ''''') ;
    end
    if ismatrix(im) ,
        % Display as a grayscale image
        if isa(im, 'uint8') && (isempty(scale) || strcmp(scale, 'dont')) ,
            % Display with no scaling
            f = figure() ;
            a = axes(f) ;
            ih = image(a, im) ;
            ih.CDataMapping = 'direct' ;
            axis(a, 'equal') ;
            axis(a, 'tight') ;            
            %clim(a, [0 255])
            colormap(a, 'gray') ;
            %colorbar(a) ;
        else
            lo = min(min(im)) ;
            hi = max(max(im)) ;
            im = double(im) ;
            f = figure() ;
            a = axes(f) ;
            ih = image(a, im) ;
            ih.CDataMapping = 'scaled' ;
            axis(a, 'equal') ;
            axis(a, 'tight') ;            
            if isempty(scale) ,
                if 0<=lo && hi<=1 ,
                    clim(a, [0 1]) ;
                    do_add_colorbar = true ;
                else
                    clim(a, [lo hi]) ;
                    do_add_colorbar = true ;
                end
            elseif strcmp(scale, 'scale')
                clim(a, [lo hi]) ;
                do_add_colorbar = true ;
            elseif strcmp(scale, 'dont')
                clim(a, [0 1]) ;
                do_add_colorbar = true ;
            else
                error('Internal error: scale has an illegal value, ''%s''', scale) ;
            end
            colormap(a, 'gray') ;
            if do_add_colorbar ,
                colorbar(a) ;
            end
        end
    else
        if ndims(im) > 3 ,
            error('imglance() doesn''t work on arrays of dimension %d > 3', ndims(im)) ;
        end
        [r,c,p] = size(im) ;
        if p == 1 ,
            % Convert to a 2D image and recurse
            im  = reshape(im, [r c]) ;
            imglance(im, scale) ;
        elseif p == 3 ,
            if isa(im, 'uint8') && (isempty(scale) || strcmp(scale, 'dont')) ,
                % Display with no scaling
                f = figure() ;
                a = axes(f) ;
                ih = image(a, im) ;
                axis(a, 'equal') ;
                axis(a, 'tight') ;
                colormap(a, 'gray') ;                
            else
                im = double(im) ;
                lo = min(min(min(im))) ;
                hi = max(max(max(im))) ;
                if isempty(scale)  ,
                    if 0<=lo && hi<=1 ,
                        % Don't scale pels
                    else
                        % Check if the pixel values are ints on [0,255], otherwise scale.
                        if does_look_like_uint8(im, lo, hi) ,                                                       
                            % If all pels are integers on [0,255], convert to uint8
                            im = uint8(im) ;
                        else
                            % Scale all pixels, channels to be on [0,1]
                            % Scale using the min, max values in the image
                            im = (im-lo)/(hi-lo) ;
                        end
                    end
                elseif strcmp(scale, 'scale')
                    % Scale all pixels, channels to be on [0,1]
                    im = (im-lo)/(hi-lo) ;
                elseif strcmp(scale, 'dont')
                    % Don't scale pixels
                else
                    error('Internal error: scale has an illegal value, ''%s''', scale) ;
                end
                f = figure() ;
                a = axes(f) ;
                ih = image(a, im) ;
                axis(a, 'equal') ;
                axis(a, 'tight') ;
            end
        else
            error('imglance() doesn''t work on 3D arrays with %d pages.  Page count must be 1 or 3.', p) ;
        end
    end
end



function result = does_look_like_uint8(im, lo, hi)
    % Test if a double mxnx3 image has all integer pixel values on [0,255]
    if 0<=lo && hi<=255 ,
        result = all(all(all(im==round(im)))) ;
    else
        result = false ;
    end
end
