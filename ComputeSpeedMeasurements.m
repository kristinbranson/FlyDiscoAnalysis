% obj.ComputeSpeedMeasurements_trx(trx,[filename])

function [trx,units] = ComputeSpeedMeasurements(trx,params,speed_file)

nflies = length(trx);
fns = {'dtheta','du_ctr','dv_ctr',...
  'corisonfly','corfrac_maj','corfrac_min','abscorfrac_min','du_cor',...
  'dv_cor','velmag_ctr','velmag','accmag','signdtheta','absdv_cor',...
  'flipdv_cor','absdtheta','d2theta','absd2theta','smooththeta',...
  'smoothdtheta','abssmoothdtheta','smoothd2theta','abssmoothd2theta',...
  'phi','yaw','absyaw',...
  'du_tail','dv_tail','absdu_tail','absdv_tail','dtheta_tail',...
  'absdtheta_tail','phisideways','absphisideways'};

isfile = exist('speed_file','var');

if isfile && exist(speed_file,'file'),
  fprintf('Loading speed measurements from file %s\n',speed_file);
  load(speed_file,'speed','units');
else
  
  fprintf('Computing speed measurements for file %s\n',speed_file);
  speed = structallocate(fns,[1,nflies]);
  units = struct;
  
  for fly = 1:nflies,
    
    % change in body orientation
    speed(fly).dtheta = modrange(diff(trx(fly).theta_mm),-pi,pi)./trx(fly).dt;
    units.dtheta = parseunits('rad/s');
    
    % change in center position
    dx = diff(trx(fly).x_mm);
    dy = diff(trx(fly).y_mm);
    
    % forward motion of body center
    if trx(fly).nframes < 2,
      speed(fly).du_ctr = [];
    else
      speed(fly).du_ctr = (dx.*cos(trx(fly).theta_mm(1:end-1)) + dy.*sin(trx(fly).theta_mm(1:end-1)))./trx(fly).dt;
    end
    units.du_ctr = parseunits('mm/s');
    % sideways motion of body center
    if trx(fly).nframes < 2,
      speed(fly).dv_ctr = [];
    else
      speed(fly).dv_ctr = (dx.*cos(trx(fly).theta_mm(1:end-1)+pi/2) + dy.*sin(trx(fly).theta_mm(1:end-1)+pi/2))./trx(fly).dt;
    end
    units.dv_ctr = parseunits('mm/s');
    
    % find the center of rotation
    if trx(fly).nframes < 2,
      corfrac = zeros(2,0);
      speed(fly).corisonfly = [];
    else
      [corfrac,speed(fly).corisonfly] = center_of_rotation2(trx(fly),false);
    end
    speed(fly).corfrac_maj = corfrac(1,:);
    speed(fly).corfrac_min = corfrac(2,:);
    speed(fly).abscorfrac_min = abs(corfrac(2,:));
    units.corfrac_maj = parseunits('unit');
    units.corfrac_min = parseunits('unit');
    units.abscorfrac_min = parseunits('unit');
    units.corisonfly = parseunits('unit');
    
    if trx(fly).nframes < 2,
      speed(fly).du_cor = [];
      speed(fly).dv_cor = [];
    else
      
      [x_cor_curr,y_cor_curr,x_cor_next,y_cor_next] = rfrac2center(trx(fly),[speed(fly).corfrac_maj;speed(fly).corfrac_min]);
      
      % change in center of rotation
      dx_cor = x_cor_next - x_cor_curr;
      dy_cor = y_cor_next - y_cor_curr;
      
      % forward motion of center of rotation
      speed(fly).du_cor = (dx_cor.*cos(trx(fly).theta_mm(1:end-1)) + dy_cor.*sin(trx(fly).theta_mm(1:end-1)))./trx(fly).dt;
      % sideways motion of body center
      speed(fly).dv_cor = (dx_cor.*cos(trx(fly).theta_mm(1:end-1)+pi/2) + dy_cor.*sin(trx(fly).theta_mm(1:end-1)+pi/2))./trx(fly).dt;
    end
    units.du_cor = parseunits('mm/s');
    units.dv_cor = parseunits('mm/s');
    
    % magnitude of velocity
    if trx(fly).nframes < 2,
      speed(fly).velmag_ctr = [];
    else
      speed(fly).velmag_ctr = sqrt(dx.^2 + dy.^2)./trx(fly).dt;
    end
    units.velmag_ctr = parseunits('mm/s');
    if trx(fly).nframes < 2,
      speed(fly).velmag = [];
    else
      speed(fly).velmag = sqrt(dx_cor.^2 + dy_cor.^2)./trx(fly).dt;
      badidx = isnan(dx_cor);
      speed(fly).velmag(badidx) = speed(fly).velmag_ctr(badidx);
    end
    units.velmag = parseunits('mm/s');
    
    % acceleration magnitude
    if trx(fly).nframes < 2,
      speed(fly).accmag = [];
    elseif trx(fly).nframes == 2,
      speed(fly).accmag = 0;
    else
      % speed from frame 1 to 2 minus speed from 2 to 3 / time from 2 to 3
      tmp = sqrt(diff(dx./trx(fly).dt).^2 + diff(dy./trx(fly).dt).^2)./trx(fly).dt(2:end);
      speed(fly).accmag = [0,tmp];
    end
    units.accmag = parseunits('mm/s/s');
    
    % flipped sign dv, dtheta
    speed(fly).signdtheta = sign(speed(fly).dtheta);
    units.signdtheta = parseunits('unit');
    speed(fly).absdv_cor = abs(speed(fly).dv_cor);
    units.absdv_cor = parseunits('mm/s');
    speed(fly).flipdv_cor = speed(fly).dv_cor.*speed(fly).signdtheta;
    units.flipdv_cor = parseunits('mm/s');
    %speed(fly).realabsdv_cor = abs(speed(fly).dv_cor);
    speed(fly).absdtheta = abs(speed(fly).dtheta);
    units.absdtheta = parseunits('rad/s');
    if trx(fly).nframes < 2,
      speed(fly).d2theta = [];
    elseif trx(fly).nframes == 2,
      speed(fly).d2theta = 0;
    else
      speed(fly).d2theta = [0,modrange(diff(speed(fly).dtheta),-pi,pi)]./trx(fly).dt;
    end
    units.d2theta = parseunits('rad/s/s');
    speed(fly).absd2theta = abs(speed(fly).d2theta);
    units.absd2theta = parseunits('rad/s/s');
    
    % smoothed orientation
    speed(fly).smooththeta = myconv(unwrap(trx(fly).theta_mm),params.thetafil,'replicate','same');
    units.smooththeta = parseunits('rad');
    if trx(fly).nframes < 2,
      speed(fly).smoothdtheta = [];
    else
      speed(fly).smoothdtheta = diff(speed(fly).smooththeta)./trx(fly).dt;
    end
    units.smoothdtheta = parseunits('rad/s');
    speed(fly).smooththeta = modrange(speed(fly).smooththeta,-pi,pi);
    speed(fly).abssmoothdtheta = abs(speed(fly).smoothdtheta);
    units.abssmoothdtheta = parseunits('rad/s');
    if trx(fly).nframes < 2,
      speed(fly).smoothd2theta = [];
    elseif trx(fly).nframes == 2,
      speed(fly).smoothd2theta = 0;

    else
      speed(fly).smoothd2theta = [0,modrange(diff(speed(fly).smoothdtheta),-pi,pi)]./trx(fly).dt;
    end
    units.smoothd2theta = parseunits('rad/s/s');
    speed(fly).abssmoothd2theta = abs(speed(fly).smoothd2theta);
    units.abssmoothd2theta = parseunits('rad/s/s');
    
    %trx(fly).off = -trx(fly).firstframe + 1;
    %trx(fly).f2i = @(f) f - trx(fly).firstframe + 1;
    
    % velocity direction
    if trx(fly).nframes < 2,
      % if only one frame, set to orientation
      speed(fly).phi = trx(fly).theta_mm;
    else
      dy1 = [trx(fly).y(2)-trx(fly).y(1),(trx(fly).y(3:end)-trx(fly).y(1:end-2))/2,trx(fly).y(end)-trx(fly).y(end-1)];
      dx1 = [trx(fly).x(2)-trx(fly).x(1),(trx(fly).x(3:end)-trx(fly).x(1:end-2))/2,trx(fly).x(end)-trx(fly).x(end-1)];
      speed(fly).phi = atan2(dy1,dx1);
    end
    units.phi = parseunits('rad');
    
    % difference between velocity direction and orientation
    speed(fly).yaw = modrange(speed(fly).phi - trx(fly).theta_mm,-pi,pi);
    units.yaw = parseunits('rad');
    speed(fly).absyaw = abs(speed(fly).yaw);
    units.absyaw = parseunits('rad');
    
    % crabwalk features
    
    % location of tail
    tailx = trx(fly).x_mm + 2*cos(-trx(fly).theta).*trx(fly).a_mm;
    taily = trx(fly).y_mm + 2*sin(-trx(fly).theta).*trx(fly).a_mm;
    % location of nose
    nosex = trx(fly).x_mm + 2*cos(trx(fly).theta).*trx(fly).a_mm;
    nosey = trx(fly).y_mm + 2*sin(trx(fly).theta).*trx(fly).a_mm;

    dx = diff(tailx);
    dy = diff(taily);

    % project onto body coords
    if trx(fly).nframes < 2,
      trx(fly).du_tail = [];
      trx(fly).dv_tail = [];
    else
      trx(fly).du_tail = dx.*cos(trx(fly).theta(1:end-1)) + dy.*sin(trx(fly).theta(1:end-1))./trx(fly).dt;
      trx(fly).dv_tail = dx.*cos(trx(fly).theta(1:end-1)+pi/2) + dy.*sin(trx(fly).theta(1:end-1)+pi/2)./trx(fly).dt;
    end
    trx(fly).units.du_tail = parseunits('mm/s');
    trx(fly).units.dv_tail = parseunits('mm/s');
    trx(fly).absdu_tail = abs(trx(fly).du_tail);
    trx(fly).units.absdu_tail = parseunits('mm/s');
    trx(fly).absdv_tail = abs(trx(fly).dv_tail);
    trx(fly).units.absdv_tail = parseunits('mm/s');
    
    % compute the rotation of nose around mean tail location

    if trx(fly).nframes < 2,
      trx(fly).dtheta_tail = [];
    else
      meantailx = (tailx(1:end-1)+tailx(2:end))/2;
      meantaily = (taily(1:end-1)+taily(2:end))/2;
      anglenose1 = atan2(nosey(1:end-1)-meantaily,nosex(1:end-1)-meantailx);
      anglenose2 = atan2(nosey(2:end)-meantaily,nosex(2:end)-meantailx);
      trx(fly).dtheta_tail = modrange(anglenose2-anglenose1,-pi,pi)./trx(fly).dt;
    end
    trx(fly).units.dtheta_tail = parseunits('rad/s');
    trx(fly).absdtheta_tail = abs(trx(fly).dtheta_tail);
    trx(fly).units.absdtheta_tail = parseunits('rad/s');
    
    % how sideways is the velocity direction?
    if trx(fly).nframes < 2,
      trx(fly).phisideways = [];
    else
      phi = atan2(dy,dx);
      trx(fly).phisideways = modrange(phi-trx(fly).theta(1:end-1),-pi/2,pi/2);
    end
    trx(fly).units.phisideways = parseunits('rad');
    trx(fly).absphisideways = abs(trx(fly).phisideways);
    trx(fly).units.absphisideways = parseunits('rad');
    
  end
  
  
  
  
  if isfile,
    fprintf('Saving speed measurements to file %s\n',speed_file);
    save(speed_file,'speed','units');
  end
end

% add these fields to trx
for i = 1:length(fns),
  fn = fns{i};
  for fly = 1:nflies,
    trx(fly).(fn) = speed(fly).(fn);
  end
end
