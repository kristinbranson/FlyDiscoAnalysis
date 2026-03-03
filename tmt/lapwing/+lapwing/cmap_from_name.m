function cmap = cmap_from_name(cmap_name)
    % set the colormap
    if strcmp(cmap_name,'red_green') ,
        % feval doesn't work with imported functions?
        cmap=lapwing.red_green(256);
    elseif strcmp(cmap_name,'red_blue') ,
        % feval doesn't work with imported functions?
        cmap=lapwing.red_blue(256);
    elseif strcmp(cmap_name,'parula') && verLessThan('matlab','8.4') ,
        cmap = jet(256) ;  % no parula colormap in early versions
    else
        cmap=feval(cmap_name,256);
    end
end    
