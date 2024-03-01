function trx = appendPhysicalUnitFieldsToTrx(input_trx, registration_data, are_timestamps_reliable, fallback_dt)
    % Add the x_mm, y_mm, a_mm, b_mm, theta_mm, dt, fps, and pxpermm fields to the
    % input trx.
    % This is a pure function.

    registerfn = registration_data.registerfn ;
    scale = registration_data.scale ;
    trx = input_trx ;
    for fly = 1:length(trx),
        % apply transformation to center, 4 extremal points on the ellipse
        x0 = trx(fly).x ;
        y0 = trx(fly).y ;
        a0 = trx(fly).a ;
        b0 = trx(fly).b ;
        theta0 = trx(fly).theta ;
        xnose0 = x0 + 2*a0.*cos(theta0);
        ynose0 = y0 + 2*a0.*sin(theta0);
        xtail0 = x0 - 2*a0.*cos(theta0);
        ytail0 = y0 - 2*a0.*sin(theta0);
        xleft0 = x0 + 2*b0.*cos(theta0-pi/2);
        yleft0 = y0 + 2*b0.*sin(theta0-pi/2);
        xright0 = x0 + 2*b0.*cos(theta0+pi/2);
        yright0 = y0 + 2*b0.*sin(theta0+pi/2);
        [~,~,chamber_index] = registerfn(x0,y0);
        [xnose1,ynose1] = registerfn(xnose0,ynose0);
        [xtail1,ytail1] = registerfn(xtail0,ytail0);
        [xleft1,yleft1] = registerfn(xleft0,yleft0);
        [xright1,yright1] = registerfn(xright0,yright0);
        % compute the center as the mean of these four points
        x1 = (xnose1+xtail1+xleft1+xright1)/4;
        y1 = (ynose1+ytail1+yleft1+yright1)/4;
        % compute the major axis from the nose to tail distance
        a1 = sqrt( (xnose1-xtail1).^2 + (ynose1-ytail1).^2 ) / 4;
        % compute the minor axis length as the left to right distance
        b1 = sqrt( (xleft1-xright1).^2 + (yleft1-yright1).^2 ) / 4;
        % compute the orientation from the nose and tail points only
        theta1 = atan2(ynose1-ytail1,xnose1-xtail1);
        % store the registerd positions
        trx(fly).x_mm = x1;
        trx(fly).y_mm = y1;
        trx(fly).a_mm = a1;
        trx(fly).b_mm = b1;
        trx(fly).theta_mm = theta1;
        trx(fly).chamber_index = mode(chamber_index) ;  % only save one (scalar) value per trajectory

        % add dt
        % TODO: fix timestamps after fix errors!
        if are_timestamps_reliable ,
            if isfield(trx,'timestamps'),
                trx(fly).dt = diff(trx(fly).timestamps);
            else
                trx(fly).dt = repmat(1/trx(fly).fps,[1,trx(fly).nframes-1]);
            end
        else
            trx(fly).dt = repmat(fallback_dt,[1,trx(fly).nframes-1]);
        end
        trx(fly).fps = 1/mean(trx(fly).dt);
        trx(fly).pxpermm = 1 / scale ;
    end
end  % function
