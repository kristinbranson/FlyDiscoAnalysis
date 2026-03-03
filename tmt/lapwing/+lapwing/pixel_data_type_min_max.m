function [d_min, d_max] = pixel_data_type_min_max(data_type)
    if isequal(data_type, 'uint8') ,
        d_min=0;
        d_max=255;
    elseif isequal(data_type, 'uint16') ,
        d_min=0;
        d_max=65535;
    else
        % This is to handle float, double
        d_min=0;
        d_max=1;
    end            
end
