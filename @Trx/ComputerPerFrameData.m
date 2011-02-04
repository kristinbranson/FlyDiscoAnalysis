function [data,units] = ComputerPerFrameData(obj,fn,n)

data = cell(1,obj.nfliespermovie(n));
nflies = obj.nfliespermovie(n);
nframes = obj.nframes(obj.exp2flies);

switch fn,
  
  % area
  case 'area',
    for fly = 1:nflies,
      a_mm = obj.GetPerFrameData('a_mm',n,fly);
      b_mm = obj.GetPerFrameData('b_mm',n,fly);
      data{fly} = (2*a_mm).*(2*b_mm)*pi;
    end
    units = parseunits('mm^2');
  
  % change in orientation
  case 'dtheta',
    for fly = 1:nflies,
      theta_mm = obj.GetPerFrameData('theta_mm',n,fly);
      dt = obj.GetPerFrameData('dt',n,fly);
      data{fly} = modrange(diff(theta_mm),-pi,pi)./dt;
    end
    units = parseunits('rad/s');

  % change in center position
  case 'du_ctr',
    for fly = 1:nflies,
      dx = diff(obj.GetPerFrameData('x_mm',n,fly));
      dy = diff(obj.GetPerFrameData('y_mm',n,fly));
      theta_mm = obj.GetPerFrameData('theta_mm',n,fly);
      dt = obj.GetPerFrameData('dt',n,fly);
      % forward motion of body center
      if nframes(fly) < 2,
        data{fly} = [];
      else
        data{fly} = (dx.*cos(theta_mm(1:end-1)) + dy.*sin(theta_mm(1:end-1)))./dt;
      end
    end
    units = parseunits('mm/s');

  % sideways motion of body center
  case 'dv_ctr',
    for fly = 1:nflies,
      dx = diff(obj.GetPerFrameData('x_mm',n,fly));
      dy = diff(obj.GetPerFrameData('y_mm',n,fly));
      theta_mm = obj.GetPerFrameData('theta_mm',n,fly);
      dt = obj.GetPerFrameData('dt',n,fly);
      % sideways motion of body center
      if nframes(fly) < 2,
        data{fly} = [];
      else
        data{fly} = (dx.*cos(theta_mm(1:end-1)+pi/2) + dy.*sin(theta_mm(1:end-1)+pi/2))./dt;
      end
    end
    units = parseunits('mm/s');

  % center of rotation
  case 'corfrac_maj',
    for fly = 1:nflies,
      if nframes(fly) < 2,
        data{fly} = [];
      else
        corfrac = obj.center_of_rotation2(n,fly,false);
        data{fly} = corfrac(1,:);
      end
    end
    units = parseunits('unit');
    
  case 'corfrac_min',
    for fly = 1:nflies,
      if nframes(fly) < 2,
        data{fly} = [];
      else
        corfrac = obj.center_of_rotation2(n,fly,false);
        data{fly} = corfrac(2,:);
      end
    end
    units = parseunits('unit');
    
  case 'corisonfly',
    for fly = 1:nflies,
      if nframes(fly) < 2,
        data{fly} = [];
      else
        [~,data{fly}] = obj.center_of_rotation2(n,fly,false);
      end
    end
    units = parseunits('unit');
    
  case 'abscorfrac_min',
    for fly = 1:nflies,
      data{fly} = abs(obj.GetPerFrameData('corfrac_min',n,fly));
    end
    units = parseunits('unit');

  % forward motion of center of rotation
  case 'du_cor',
    for fly = 1:nflies,
      if nframes(fly) < 2,
        data{fly} = [];
      else
        [x_cor_curr,y_cor_curr,x_cor_next,y_cor_next] = ...
          obj.rfrac2center(n,fly);
        theta_mm = obj.GetPerFrameData('theta_mm',n,fly);
        dt = obj.GetPerFrameData('dt',n,fly);
        dx_cor = x_cor_next - x_cor_curr;
        dy_cor = y_cor_next - y_cor_curr;
        data{fly} = (dx_cor.*cos(theta_mm(1:end-1)) + ...
          dy_cor.*sin(theta_mm(1:end-1)))./dt;
      end
    end
    units = parseunits('mm/s');

  % sideways motion of center of rotation
  case 'dv_cor',
    for fly = 1:nflies,
      if nframes(fly) < 2,
        data{fly} = [];
      else
        [x_cor_curr,y_cor_curr,x_cor_next,y_cor_next] = ...
          obj.rfrac2center(n,fly);
        theta_mm = obj.GetPerFrameData('theta_mm',n,fly);
        dt = obj.GetPerFrameData('dt',n,fly);
        dx_cor = x_cor_next - x_cor_curr;
        dy_cor = y_cor_next - y_cor_curr;
        data{fly} = (dx_cor.*cos(theta_mm(1:end-1)+pi/2) + ...
          dy_cor.*sin(theta_mm(1:end-1)+pi/2))./dt;
      end
    end
    units = parseunits('mm/s');

  % magnitude of velocity
  case 'velmag_ctr',
    for fly = 1:nflies,
      if nframes(fly) < 2,
        data{fly} = [];
      else
        dx = diff(obj.GetPerFrameData('x_mm',n,fly));
        dy = diff(obj.GetPerFrameData('y_mm',n,fly));
        dt = obj.GetPerFrameData('dt',n,fly);
        data{fly} = sqrt(dx.^2 + dy.^2)./dt;
      end
    end
    units = parseunits('mm/s');
    
  % magnitude of velocity around center of rotation
  case 'velmag',
    for fly = 1:nflies,
      if nframes(fly) < 2,
        data{fly} = [];
      else
        [x_cor_curr,y_cor_curr,x_cor_next,y_cor_next] = ...
          obj.rfrac2center(n,fly);
        dt = obj.GetPerFrameData('dt',n,fly);
        dx_cor = x_cor_next - x_cor_curr;
        dy_cor = y_cor_next - y_cor_curr;
        data{fly} = sqrt(dx_cor.^2 + dy_cor.^2)./dt;
        badidx = isnan(dx_cor);
        velmag_ctr = obj.GetPerFrameData('velmag_ctr',n,fly);
        data{fly}(badidx) = velmag_ctr(badidx);
      end
    end
    
    units = parseunits('mm/s');

  % acceleration magnitude
  case 'accmag',
    for fly = 1:nflies,
      if nframes(fly) < 2,
        data{fly} = [];
      elseif nframes(fly) == 2,
        data{fly} = 0;
      else
        % speed from frame 1 to 2 minus speed from 2 to 3 / time from 2 to 3
        dx = diff(obj.GetPerFrameData('x_mm',n,fly));
        dy = diff(obj.GetPerFrameData('y_mm',n,fly));
        dt = obj.GetPerFrameData('dt',n,fly);
        tmp = sqrt(diff(dx./dt).^2 + diff(dy./dt).^2)./dt(2:end);
        data{fly} = [0,tmp];
      end
    end
    units = parseunits('mm/s/s');

  % flipped sign dv, dtheta
  case 'signdtheta',
    for fly = 1:nflies,
      data{fly} = sign(obj.GetPerFrameData('dtheta',n,fly));
    end
    units = parseunits('unit');
    
  case 'absdv_cor',
    for fly = 1:nflies,
      data{fly} = abs(obj.GetPerFrameData('dv_cor',n,fly));
    end
    units = parseunits('mm/s');
    
  case 'flipdv_cor',
    for fly = 1:nflies,
      data{fly} = obj.GetPerFrameData('dv_cor',n,fly).* ...
        obj.GetPerFrameData('signdtheta',n,fly);
    end
    units = parseunits('mm/s');

  case 'absdtheta',
    for fly = 1:nflies,
      data{fly} = abs(obj.GetPerFrameData('dtheta',n,fly));
    end
    units = parseunits('rad/s');
    
  case 'd2theta',
    for fly = 1:nflies,
      if nframes(fly) < 2,
        data{fly} = [];
      elseif nframes(fly) == 2,
        data{fly} = 0;
      else
        dt = obj.GetPerFrameData('dt',n,fly);
        dtheta = obj.GetPerFrameData('dtheta',n,fly);
        data{fly} = [0,modrange(diff(dtheta),-pi,pi)]./dt;
      end
    end
    units = parseunits('rad/s/s');
    
    case 'absd2theta',
      for fly = 1:nflies,
        data{fly} = abs(obj.GetPerFrameData('d2theta',n,fly));
      end
      units = parseunits('rad/s/s');

    % smoothed orientation
  case 'smooththeta',
    for fly = 1:nflies,
      theta_mm = obj.GetPerFrameData('theta_mm',n,fly);
      data{fly} = modrange(myconv(unwrap(theta_mm),obj.thetafil,'replicate','same'),-pi,pi);
    end
    units = parseunits('rad');
    
  case 'smoothdtheta',
    for fly = 1:nflies,
      if nframes(fly) < 2,
        data{fly} = [];
      else
        theta_mm = obj.GetPerFrameData('theta_mm',n,fly);
        smooththeta = myconv(unwrap(theta_mm),obj.thetafil,'replicate','same');
        dt = obj.GetPerFrameData('dt',n,fly);
        data{fly} = diff(smooththeta)./dt;
      end
    end
    units = parseunits('rad/s');
    
  case 'abssmoothdtheta',
    for fly = 1:nflies,
      data{fly} = abs(obj.GetPerFrameData('smoothdtheta',n,fly));
    end
    units = parseunits('rad/s');

  case 'smoothd2theta',
    for fly = 1:nflies,
      if nframes(fly) < 2,
        data{fly} = [];
      elseif nframes(fly) == 2,
        data{fly} = 0;
      else
        smoothdtheta = obj.GetPerFrameData('smoothdtheta',n,fly);
        dt = obj.GetPerFrameData('dt',n,fly);
        data{fly} = [0,modrange(diff(smoothdtheta),-pi,pi)]./dt;
      end
    end
    units = parseunits('rad/s/s');
    
  case 'abssmoothd2theta',
    for fly = 1:nflies,
      data{fly} = abs(obj.GetPerFrameData,'smoothd2theta',n,fly);
    end
    units = parseunits('rad/s/s');
      
  % velocity direction
  case 'phi',
    for fly = 1:nflies,
      if nframes(fly) < 2,
        % if only one frame, set to orientation
        data{fly} = obj.GetPerFrameData('theta_mm',n,fly);
      else
        x_mm = obj.GetPerFrameData('x_mm',n,fly);
        y_mm = obj.GetPerFrameData('y_mm',n,fly);
        dy1 = [y_mm(2)-y_mm(1),(y_mm(3:end)-y_mm(1:end-2))/2,y_mm(end)-y_mm(end-1)];
        dx1 = [x_mm(2)-x_mm(1),(x_mm(3:end)-x_mm(1:end-2))/2,x_mm(end)-x_mm(end-1)];
        data{fly} = atan2(dy1,dx1);
      end
    end
    units = parseunits('rad');
    
  % difference between velocity direction and orientation
  case 'yaw',
    for fly = 1:nflies,
      phi = obj.GetPerFrameData('phi',n,fly);
      theta_mm = obj.GetPerFrameData('theta_mm',n,fly);
      data{fly} = modrange(phi - theta_mm,-pi,pi);
    end
    units = parseunits('rad');
    
  case 'absyaw',
    for fly = 1:nflies,
      data{fly} = abs(obj.GetPerFrameData('yaw',n,fly));
    end
    units = parseunits('rad');
      
    
    
  % crabwalk features
  case 'du_tail',
    
    for fly = 1:nflies,
      
      if nframes(fly) < 2,
        data{fly} = [];
      else
        x_mm = obj.GetPerFrameData('x_mm',n,fly);
        y_mm = obj.GetPerFrameData('y_mm',n,fly);
        theta_mm = obj.GetPerFrameData('theta_mm',n,fly);
        a_mm = obj.GetPerFrameData('a_mm',n,fly);
        dt = obj.GetPerFrameData('dt',n,fly);
        
        % location of tail
        tailx = x_mm + 2*cos(-theta_mm).*a_mm;
        taily = y_mm + 2*sin(-theta_mm).*a_mm;
         
        dx = diff(tailx);
        dy = diff(taily);

        data{fly} = dx.*cos(theta_mm(1:end-1)) + dy.*sin(theta_mm(1:end-1))./dt;
        
      end
    end
    units = parseunits('mm/s');
    
  case 'dv_tail',
        
    for fly = 1:nflies,
      
      if nframes(fly) < 2,
        data{fly} = [];
      else
        x_mm = obj.GetPerFrameData('x_mm',n,fly);
        y_mm = obj.GetPerFrameData('y_mm',n,fly);
        theta_mm = obj.GetPerFrameData('theta_mm',n,fly);
        a_mm = obj.GetPerFrameData('a_mm',n,fly);
        dt = obj.GetPerFrameData('dt',n,fly);
        
        % location of tail
        tailx = x_mm + 2*cos(-theta_mm).*a_mm;
        taily = y_mm + 2*sin(-theta_mm).*a_mm;
        
        dx = diff(tailx);
        dy = diff(taily);

        data{fly} = dx.*cos(theta_mm(1:end-1)+pi/2) + dy.*sin(theta_mm(1:end-1)+pi/2)./dt;
        
      end
    end
    units = parseunits('mm/s');

  case 'absdu_tail',
    for fly = 1:nflies,
      data{fly} = abs(obj.GetPerFrameData('du_tail',n,fly));
    end
    units = parseunits('mm/s');

  case 'absdv_tail',
    for fly = 1:nflies,
      data{fly} = abs(obj.GetPerFrameData('dv_tail',n,fly));
    end
    units = parseunits('mm/s');

  % compute the rotation of nose around mean tail location
  case 'dtheta_tail',
    for fly = 1:nflies,
      if nframes(fly) < 2,
        data{fly} = [];
      else

        x_mm = obj.GetPerFrameData('x_mm',n,fly);
        y_mm = obj.GetPerFrameData('y_mm',n,fly);
        theta_mm = obj.GetPerFrameData('theta_mm',n,fly);
        a_mm = obj.GetPerFrameData('a_mm',n,fly);
        dt = obj.GetPerFrameData('dt',n,fly);
        
        % location of tail
        tailx = x_mm + 2*cos(-theta_mm).*a_mm;
        taily = y_mm + 2*sin(-theta_mm).*a_mm;
        
        % location of nose
        nosex = x_mm + 2*cos(theta_mm).*a_mm;
        nosey = y_mm + 2*sin(theta_mm).*a_mm;

        meantailx = (tailx(1:end-1)+tailx(2:end))/2;
        meantaily = (taily(1:end-1)+taily(2:end))/2;
        anglenose1 = atan2(nosey(1:end-1)-meantaily,nosex(1:end-1)-meantailx);
        anglenose2 = atan2(nosey(2:end)-meantaily,nosex(2:end)-meantailx);
        data{fly} = modrange(anglenose2-anglenose1,-pi,pi)./dt;
      end
    end
    units = parseunits('rad/s');
    
  case 'absdtheta_tail',
    for fly = 1:nflies,
      data{fly} = abs(obj.GetPerFrameData('dtheta_tail',n,fly));
    end
    units = parseunits('rad/s');

  % how sideways is the velocity direction?
  case 'phisideways',
    for fly = 1:nflies,
      if nframes(fly) < 2,
        data{fly} = [];
      else
        x_mm = obj.GetPerFrameData('x_mm',n,fly);
        y_mm = obj.GetPerFrameData('y_mm',n,fly);
        theta_mm = obj.GetPerFrameData('theta_mm',n,fly);
        a_mm = obj.GetPerFrameData('a_mm',n,fly);
         
        % location of tail
        tailx = x_mm + 2*cos(-theta_mm).*a_mm;
        taily = y_mm + 2*sin(-theta_mm).*a_mm;
        dx = diff(tailx);
        dy = diff(taily);

        phi = atan2(dy,dx);
        data{fly} = modrange(phi-theta_mm(1:end-1),-pi/2,pi/2);
      end
    end
    units = parseunits('rad');
    
  case 'absphisideways',
    for fly = 1:nflies,
      data{fly} = abs(obj.GetPerFrameData('phisideways',n,fly));
    end
    units = parseunits('rad');

  % distance from center
  case 'arena_r',
    for fly = 1:nflies,
      x_mm = obj.GetPerFrameData('x_mm',n,fly);
      y_mm = obj.GetPerFrameData('y_mm',n,fly);
      wall = obj.GetArenaWallParams(n);
      data{fly} = sqrt((x_mm - wall.arena_center_mm(1)).^2 + ...
        (y_mm - wall.arena_center_mm(2)).^2);
    end
    units = parseunits('mm');
    
  % distance to wall
  case 'dist2wall',
    for fly = 1:nflies,
      arena_r = obj.GetPerFrameData('arena_r',n,fly);
      wall = obj.GetArenaWallParams(n);
      data{fly} = wall.arena_radius_mm - arena_r;
    end
    units = parseunits('mm');
    
  % polar angle of closest point on the wall
  case 'wallangle',

    for fly = 1:nflies,
      x_mm = obj.GetPerFrameData('x_mm',n,fly);
      y_mm = obj.GetPerFrameData('y_mm',n,fly);
      wall = obj.GetArenaWallParams(n);
      data{fly} = atan2(y_mm-wall.arena_center_mm(2),...
        x_mm-wall.arena_center_mm(1));
    end
    units = parseunits('rad');
    
  % angle to closest point on the wall in the fly's coordinate system
  case 'angle2wall',
    for fly = 1:nflies,
      wallangle = obj.GetPerFrameData('wallangle',n,fly);
      theta_mm = obj.GetPerFrameData('theta_mm',n,fly);
      data{fly} = modrange(wallangle-theta_mm,-pi,pi);
    end
    units = parseunits('rad');
    
  % change in distance to wall
  case 'ddist2wall',
    for fly = 1:nflies,
      dist2wall = obj.GetPerFrameData('dist2wall',n,fly);
      dt = obj.GetPerFrameData('dt',n,fly);
      data{fly} = diff(dist2wall)./dt;
    end
    units = parseunits('mm/s');

  % change in polar angle to closest point on wall
  case 'dwallangle',
    for fly = 1:nflies,
      wallangle = obj.GetPerFrameData('wallangle',n,fly);
      dt = obj.GetPerFrameData('dt',n,fly);
      data{fly} = modrange(diff(wallangle),-pi,pi)./dt;
    end
    units = parseunits('rad/s');

  % absolute change in polar angle to closest point on wall
  case 'absdwallangle',
    for fly = 1:nflies,
      dwallangle = obj.GetPerFrameData('dwallangle',n,fly);
      data{fly} = abs(dwallangle);
    end
    units = parseunits('rad/s');
    
  % change in angle to closest point on wall in fly's coordinate system
  case 'dangle2wall',
    for fly = 1:nflies,
      angle2wall = obj.GetPerFrameData('angle2wall',n,fly);
      dt = obj.GetPerFrameData('dt',n,fly);
      data{fly} = modrange(diff(angle2wall),-pi,pi)./dt;
    end
    units = parseunits('rad/s');
    
  % position of nose
  case 'xnose_mm',

    for fly = 1:nflies,
      x_mm = obj.GetPerFrameData('x_mm',n,fly);
      a_mm = obj.GetPerFrameData('a_mm',n,fly);
      theta_mm = obj.GetPerFrameData('theta_mm',n,fly);

      data{fly} = x_mm + 2*a_mm.*cos(theta_mm);
    end
    
    units = parseunits('mm');

  case 'ynose_mm',
    for fly = 1:nflies,
      y_mm = obj.GetPerFrameData('y_mm',n,fly);
      a_mm = obj.GetPerFrameData('a_mm',n,fly);
      theta_mm = obj.GetPerFrameData('theta_mm',n,fly);

      data{fly} = y_mm + 2*a_mm.*sin(theta_mm);
    end    
    units = parseunits('mm');

  case 'closestfly_center',
    
    mindcenter = cell(1,nflies);
    closestfly = cell(1,nflies);
    for fly = 1:nflies,
      
      dcenter = obj.ComputeDistFlyCenter(n,fly);
      [mindcenter{fly},closestfly{fly}] = min(dcenter,[],1);
    end

    % so that we don't compute dcenter twice
    data = dcenter; %#ok<NASGU>
    units = parseunits('mm'); %#ok<NASGU>
    filename = obj.GetPerFrameFile('dcenter',n);
    save(filename,'data','units');

    data = closestfly;
    units = parseunits('unit');
    
  case 'dcenter',
      
    mindcenter = cell(1,nflies);
    closestfly = cell(1,nflies);
    for fly = 1:nflies,
      
      dcenter = obj.ComputeDistFlyCenter(n,fly);
      [mindcenter{fly},closestfly{fly}] = min(dcenter,[],1);

    end

    % so that we don't compute dcenter twice
    data = closestfly; %#ok<NASGU>
    units = parseunits('units'); %#ok<NASGU>
    filename = obj.GetPerFrameFile('closestfly_center',n);
    save(filename,'data','units');

    data = mindcenter;
    units = parseunits('mm');
    
    
  case 'closestfly_nose2ell',
    
    mind = cell(1,nflies);
    closestfly = cell(1,nflies);
    for fly = 1:nflies,
      
      d = obj.ComputeDistFlyNose2Ell(n,fly);
      [mind{fly},closestfly{fly}] = min(d,[],1);
    end

    % so that we don't compute dcenter twice
    data = d; %#ok<NASGU>
    units = parseunits('mm'); %#ok<NASGU>
    filename = obj.GetPerFrameFile('dnose2ell',n);
    save(filename,'data','units');

    data = closestfly;
    units = parseunits('unit');

  case 'dnose2ell',
    
    mind = cell(1,nflies);
    closestfly = cell(1,nflies);
    for fly = 1:nflies,      
      d = obj.ComputeDistFlyNose2Ell(n,fly);
      [mind{fly},closestfly{fly}] = min(d,[],1);
    end

    % so that we don't compute dcenter twice
    data = closestfly; %#ok<NASGU>
    units = parseunits('units'); %#ok<NASGU>
    filename = obj.GetPerFrameFile('closestfly_nose2ell',n);
    save(filename,'data','units');

    data = mind;
    units = parseunits('mm');
    
  case 'closestfly_ell2nose',
    
    mind = cell(1,nflies);
    closestfly = cell(1,nflies);
    for fly = 1:nflies,
      
      d = obj.ComputeDistFlyEll2Nose(n,fly);
      [mind{fly},closestfly{fly}] = min(d,[],1);
    end

    % so that we don't compute dcenter twice
    data = d; %#ok<NASGU>
    units = parseunits('mm'); %#ok<NASGU>
    filename = obj.GetPerFrameFile('dell2nose',n);
    save(filename,'data','units');

    data = closestfly;
    units = parseunits('unit');

  case 'dell2nose',
    
    mind = cell(1,nflies);
    closestfly = cell(1,nflies);
    for fly = 1:nflies,      
      d = obj.ComputeDistFlyEll2Nose(n,fly);
      [mind{fly},closestfly{fly}] = min(d,[],1);
    end

    % so that we don't compute dcenter twice
    data = closestfly; %#ok<NASGU>
    units = parseunits('units'); %#ok<NASGU>
    filename = obj.GetPerFrameFile('closestfly_ell2nose',n);
    save(filename,'data','units');

    data = mind;
    units = parseunits('mm');
    
    
  case 'closestfly_anglesub',
    
    maxanglesub = cell(1,nflies);
    closestfly = cell(1,nflies);
    for fly = 1:nflies,
      
      anglesub = obj.ComputeFlyAngleSubtended(n,fly);
      [maxanglesub{fly},closestfly{fly}] = max(anglesub,[],1);
      
    end

    % so that we don't compute dcenter twice
    data = maxanglesub; %#ok<NASGU>
    units = parseunits('rad'); %#ok<NASGU>
    filename = obj.GetPerFrameFile('anglesub',n);
    save(filename,'data','units');

    data = closestfly;
    units = parseunits('unit');

  case 'anglesub',
    
    maxanglesub = cell(1,nflies);
    closestfly = cell(1,nflies);
    for fly = 1:nflies,
      
      anglesub = obj.ComputeFlyAngleSubtended(n,fly);
      [maxanglesub{fly},closestfly{fly}] = max(anglesub,[],1);
      
    end

    % so that we don't compute dcenter twice
    
    data = closestfly; %#ok<NASGU>
    units = parseunits('unit'); %#ok<NASGU>
    filename = obj.GetPerFrameFile('anglesub',n);
    save(filename,'data','units');

    data = maxanglesub; 
    units = parseunits('rad');
    
  case {'magveldiff_center','magveldiff_nose2ell',...
      'magveldiff_ell2nose','magveldiff_anglesub'},

    for fly1 = 1:nflies,
      switch fn,
        case 'magveldiff_center',
          closestfly = obj.GetPerFrameData('closestfly_center',n,fly1);
        case 'magveldiff_nose2ell',
          closestfly = obj.GetPerFrameData('closestfly_nose2ell',n,fly1);
        case 'magveldiff_ell2nose',
          closestfly = obj.GetPerFrameData('closestfly_ell2nose',n,fly1);
        case 'magveldiff_anglesub',
          closestfly = obj.GetPerFrameData('closestfly_anglesub',n,fly1);
      end
      x_mm1 = obj.GetPerFrameData('x_mm',n,fly1);
      y_mm1 = obj.GetPerFrameData('y_mm',n,fly1);
      dt = obj.GetPerFrameData('dt',n,fly1);
      flyidx1 = obj.getFlyIdx(n,fly1);
      firstframe1 = obj.firstframes(flyidx1);
      data{fly1} = nan(1,nframes(fly1));
      
      for fly2 = 1:nflies,
        if fly1 == fly2, continue; end
        idx = find(closestfly == fly2);
        if isempty(idx), continue; end
        flyidx2 = obj.getFlyIdx(n,fly2);
        firstframe2 = obj.firstframes(flyidx2);
        off = firstframe1 - firstframe2;
        x_mm2 = obj.GetPerFrameData('x_mm',n,fly2);
        y_mm2 = obj.GetPerFrameData('y_mm',n,fly2);
        dx = x_mm2(off+idx) - x_mm1(idx);
        dy = y_mm2(off+idx) - y_mm1(idx);
        data{fly1}(idx) = sqrt(dx.^2 + dy.^2);
      end
      
      data{fly1} = data{fly1} ./ dt;
      
      units = parseunits('mm/s');
      
    end
    
  case {'veltoward_center','veltoward_nose2ell',...
      'veltoward_ell2nose','veltoward_anglesub'},

    for fly1 = 1:nflies,
      switch fn,
        case 'veltoward_center',
          closestfly = obj.GetPerFrameData('closestfly_center',n,fly1);
        case 'veltoward_nose2ell',
          closestfly = obj.GetPerFrameData('closestfly_nose2ell',n,fly1);
        case 'veltoward_ell2nose',
          closestfly = obj.GetPerFrameData('closestfly_ell2nose',n,fly1);
        case 'veltoward_anglesub',
          closestfly = obj.GetPerFrameData('closestfly_anglesub',n,fly1);
      end
      x_mm1 = obj.GetPerFrameData('x_mm',n,fly1);
      y_mm1 = obj.GetPerFrameData('y_mm',n,fly1);
      dt = obj.GetPerFrameData('dt',n,fly1);
      flyidx1 = obj.getFlyIdx(n,fly1);
      firstframe1 = obj.firstframes(flyidx1);
      data{fly1} = nan(1,nframes(fly1)-1);
      vel_x = diff(x_mm1);
      vel_y = diff(y_mm1);
      
      for fly2 = 1:nflies,
        if fly1 == fly2, continue; end
        idx = find(closestfly(1:end-1) == fly2);
        if isempty(idx), continue; end
        flyidx2 = obj.getFlyIdx(n,fly2);
        firstframe2 = obj.firstframes(flyidx2);
        off = firstframe1 - firstframe2;
        x_mm2 = obj.GetPerFrameData('x_mm',n,fly2);
        y_mm2 = obj.GetPerFrameData('y_mm',n,fly2);
        dx = x_mm2(off+idx) - x_mm1(idx);
        dy = y_mm2(off+idx) - y_mm1(idx);
        z = sqrt(dx.^2 + dy.^2);
        % direction to other fly
        dx = dx ./ z;
        dy = dy ./ z;
        
        data{fly1}(idx) = (dx.*vel_x(idx) + dy.*vel_y(idx));
      end
      
      data{fly1} = data{fly1} ./ dt;
      
      units = parseunits('mm/s');
      
    end
    
  case {'absthetadiff_center','absthetadiff_nose2ell',...
      'absthetadiff_ell2nose','absthetadiff_anglesub'},

    for fly1 = 1:nflies,
      switch fn,
        case 'absthetadiff_center',
          closestfly = obj.GetPerFrameData('closestfly_center',n,fly1);
        case 'absthetadiff_nose2ell',
          closestfly = obj.GetPerFrameData('closestfly_nose2ell',n,fly1);
        case 'absthetadiff_ell2nose',
          closestfly = obj.GetPerFrameData('closestfly_ell2nose',n,fly1);
        case 'absthetadiff_anglesub',
          closestfly = obj.GetPerFrameData('closestfly_anglesub',n,fly1);
      end
      theta_mm1 = obj.GetPerFrameData('theta_mm',n,fly1);
      flyidx1 = obj.getFlyIdx(n,fly1);
      firstframe1 = obj.firstframes(flyidx1);
      data{fly1} = nan(1,nframes(fly1)-1);
      
      for fly2 = 1:nflies,
        if fly1 == fly2, continue; end
        idx = find(closestfly == fly2);
        if isempty(idx), continue; end
        flyidx2 = obj.getFlyIdx(n,fly2);
        firstframe2 = obj.firstframes(flyidx2);
        off = firstframe1 - firstframe2;
        theta_mm2 = obj.GetPerFrameData('theta_mm',n,fly2);
        data{fly1}(idx) = abs(modrange(theta_mm2(off+idx)-theta_mm1(idx),-pi,pi));
      end
      
      units = parseunits('rad');
      
    end
    
  case {'absphidiff_center','absphidiff_nose2ell',...
      'absphidiff_ell2nose','absphidiff_anglesub'},

    for fly1 = 1:nflies,
      switch fn,
        case 'absphidiff_center',
          closestfly = obj.GetPerFrameData('closestfly_center',n,fly1);
        case 'absphidiff_nose2ell',
          closestfly = obj.GetPerFrameData('closestfly_nose2ell',n,fly1);
        case 'absphidiff_ell2nose',
          closestfly = obj.GetPerFrameData('closestfly_ell2nose',n,fly1);
        case 'absphidiff_anglesub',
          closestfly = obj.GetPerFrameData('closestfly_anglesub',n,fly1);
      end
      phi1 = obj.GetPerFrameData('phi',n,fly1);
      flyidx1 = obj.getFlyIdx(n,fly1);
      firstframe1 = obj.firstframes(flyidx1);
      data{fly1} = nan(1,nframes(fly1)-1);
      
      for fly2 = 1:nflies,
        if fly1 == fly2, continue; end
        idx = find(closestfly(1:end-1) == fly2);
        if isempty(idx), continue; end
        flyidx2 = obj.getFlyIdx(n,fly2);
        firstframe2 = obj.firstframes(flyidx2);
        off = firstframe1 - firstframe2;
        phi2 = obj.GetPerFrameData('phi',n,fly2);
        data{fly1}(idx) = abs(modrange(phi2(off+idx)-phi1(idx),-pi,pi));
      end
      
      units = parseunits('rad');
      
    end
    
  case {'absanglefrom1to2_center','absanglefrom1to2_nose2ell',...
      'absanglefrom1to2_ell2nose','absanglefrom1to2_anglesub'},
    
    for fly1 = 1:nflies,
      switch fn,
        case 'absanglefrom1to2_center',
          closestfly = obj.GetPerFrameData('closestfly_center',n,fly1);
        case 'absanglefrom1to2_nose2ell',
          closestfly = obj.GetPerFrameData('closestfly_nose2ell',n,fly1);
        case 'absanglefrom1to2_ell2nose',
          closestfly = obj.GetPerFrameData('closestfly_ell2nose',n,fly1);
        case 'absanglefrom1to2_anglesub',
          closestfly = obj.GetPerFrameData('closestfly_anglesub',n,fly1);
      end
      x_mm1 = obj.GetPerFrameData('x_mm',n,fly1);
      y_mm1 = obj.GetPerFrameData('y_mm',n,fly1);
      theta_mm1 = obj.GetPerFrameData('y_mm',n,fly1);
      flyidx1 = obj.getFlyIdx(n,fly1);
      firstframe1 = obj.firstframes(flyidx1);
      data{fly1} = nan(1,nframes(fly1)-1);
      
      for fly2 = 1:nflies,
        if fly1 == fly2, continue; end
        idx = find(closestfly == fly2);
        if isempty(idx), continue; end
        flyidx2 = obj.getFlyIdx(n,fly2);
        firstframe2 = obj.firstframes(flyidx2);
        off = firstframe1 - firstframe2;
        x_mm2 = obj.GetPerFrameData('x_mm',n,fly2);
        y_mm2 = obj.GetPerFrameData('y_mm',n,fly2);
        dx = x_mm2(off+idx) - x_mm1(idx);
        dy = y_mm2(off+idx) - y_mm1(idx);
        data{fly1}(idx) = abs(modrange(atan2(dy,dx)-theta_mm1(idx),-pi,pi));
      end
      
      units = parseunits('rad');
      
    end
    
  case {'ddcenter','ddnose2ell','ddell2nose','danglesub'},
    
    for fly = 1:nflies,
      
      switch fn,
        case 'ddcenter',
          closestfly = obj.GetPerFrameData('closestfly_center',n,fly);
          d = obj.GetPerFrameData('dcenter',n,fly);
          units = parseunits('mm/s');
        case 'ddnose2ell',
          closestfly = obj.GetPerFrameData('closestfly_nose2ell',n,fly);
          d = obj.GetPerFrameData('dnose2ell',n,fly);
          units = parseunits('mm/s');
        case 'ddell2nose',
          closestfly = obj.GetPerFrameData('closestfly_ell2nose',n,fly);
          d = obj.GetPerFrameData('dell2nose',n,fly);
          units = parseunits('mm/s');
        case 'danglesub',
          closestfly = obj.GetPerFrameData('closestfly_anglesub',n,fly);
          d = obj.GetPerFrameData('anglesub',n,fly);
          units = parseunits('rad/s');
      end
      
      issamefly = closestfly(1:end-1) == closestfly(2:end);
      dt = obj.GetPerFrameData('dt',n,fly);
      data{fly} = diff(d) ./ dt;
      data{fly}(~issamefly) = nan;
      
    end
    
  case 'sex',
    
    data = obj.ClassifySex(n);
    units = parseunits('unit');
    
  case 'areasmooth',

    f = fdesign.lowpass('N,F3db',obj.areasmooth_filterorder,obj,areasmooth_maxfreq);
    h = design(f,'butter');
    h.PersistentMemory = true;

    for fly = 1:nflies,
      area = obj.GetPerFrameData('area',n,fly);
      h.filter(fliplr(area));
      areasmooth = h.filter(area);
      isoutlier = abs(areasmooth - area) > obj.areasmooth_maxerr;
      [starts,ends] = get_interval_ends(isoutlier);
      ends = ends - 1;
      areacurr = area;
      for i = 1:numel(starts),
        if starts(i) == 1 && ends(i) == nframes(fly),
          break;
        elseif starts(i) == 1,
          areacurr(starts(i):ends(i)) = area(ends(i)+1);
        elseif ends(i) == trx(fly).nframes,
          areacurr(starts(i):ends(i)) = area(starts(i)-1);
        else
          areacurr(starts(i):ends(i)) = (area(starts(i)-1)+area(ends(i)+1))/2;
        end
      end
      data{fly} = areacurr;
    end
    
    units = parseunits('mm^2');
    
  otherwise
    error('Unknown per-frame data field name %s',fn);
end
filename = obj.GetPerFrameFile(fn,n);
save(filename,'data','units');

