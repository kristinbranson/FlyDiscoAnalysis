% obj.ComputeClosestFlyMeasurements_trx(trx,[closestfly_file])
function [trx,units] = ComputeClosestFlyMeasurements(trx,params,closestfly_file) 

nflies = length(trx);
fns = {'dcenter','closestfly_center','dnose2ell','closestfly_nose2ell',...
  'dell2nose','closestfly_ell2nose','anglesub','closestfly_anglesub',...
  'magveldiff_center','magveldiff_nose2ell','magveldiff_ell2nose',...
  'magveldiff_anglesub','veltoward_center','veltoward_nose2ell',...
  'veltoward_ell2nose','veltoward_anglesub','absthetadiff_center',...
  'absthetadiff_nose2ell','absthetadiff_ell2nose','absthetadiff_anglesub',...
  'absphidiff_center','absphidiff_nose2ell','absphidiff_ell2nose',...
  'absphidiff_anglesub','absanglefrom1to2_center',...
  'absanglefrom1to2_nose2ell','absanglefrom1to2_ell2nose',...
  'absanglefrom1to2_anglesub','ddnose2ell','ddell2nose','danglesub'};

isfile = exist('closestfly_file','var');

if isfile && exist(closestfly_file,'file'),
  fprintf('Loading closest fly measurements from file %s\n',closestfly_file);
  load(closestfly_file,'closestfly','units');
else

  fprintf('Computing closest fly measurements for file %s\n',closestfly_file);
  closestfly = structallocate(fns,[1,nflies]);
  
  for fly1 = 1:nflies,

    fprintf('First fly = %d\n',fly1);
    
    % position of nose
    xnose = trx(fly1).x_mm + 2*trx(fly1).a_mm.*cos(trx(fly1).theta_mm);
    ynose = trx(fly1).y_mm + 2*trx(fly1).a_mm.*sin(trx(fly1).theta_mm);
 
    % initialize
    dcenter = nan(nflies,trx(fly1).nframes);
    dnose2ell = nan(nflies,trx(fly1).nframes);
    dell2nose = nan(nflies,trx(fly1).nframes);
    magveldiff = nan(nflies,trx(fly1).nframes-1);
    veltoward = nan(nflies,trx(fly1).nframes-1);
    absthetadiff = nan(nflies,trx(fly1).nframes);
    absphidiff = nan(nflies,trx(fly1).nframes-1);
    absanglefrom1to2 = nan(nflies,trx(fly1).nframes);
    anglesub = nan(nflies,trx(fly1).nframes);
  
    % loop over other flies
    for fly2 = 1:nflies,
      if fly2 == fly1, continue; end
    
      fprintf('Other fly = %d\n',fly2);
    
      % get start and end frames of overlap
      t0 = max(trx(fly1).firstframe,trx(fly2).firstframe);
      t1 = min(trx(fly1).endframe,trx(fly2).endframe);
    
      % no overlap
      if t1 < t0, continue; end
    
      % indices for these frames
      offi = trx(fly1).firstframe-1;
      offj = trx(fly2).firstframe-1;
      i0 = t0 - offi;
      i1 = t1 - offi;
      j0 = t0 - offj;
      j1 = t1 - offj;

      % centroid distance
      dx = trx(fly2).x_mm(j0:j1)-trx(fly1).x_mm(i0:i1);
      dy = trx(fly2).y_mm(j0:j1)-trx(fly1).y_mm(i0:i1);
      z = sqrt(dx.^2 + dy.^2);
      dcenter(fly2,i0:i1) = z;
      % direction to other fly
      dx = dx ./ z;
      dy = dy ./ z;
    
      % distance from fly1's nose to fly2
      for t = t0:t1,
        i = t - offi;
        j = t - offj;
        dnose2ell(fly2,i) = ellipsedist_hack(trx(fly2).x_mm(j),trx(fly2).y_mm(j),...
          trx(fly2).a_mm(j),trx(fly2).b_mm(j),trx(fly2).theta_mm(j),...
          xnose(i),ynose(i));
      end
      
      % distance from fly2's nose to fly1
      for t = t0:t1,
        i = t - offi;
        j = t - offj;
        xnose2 = trx(fly2).x_mm(j) + 2*trx(fly2).a_mm(j).*cos(trx(fly2).theta_mm(j));
        ynose2 = trx(fly2).y_mm(j) + 2*trx(fly2).a_mm(j).*sin(trx(fly2).theta_mm(j));
        dell2nose(fly2,i) = ellipsedist_hack(trx(fly1).x_mm(i),trx(fly1).y_mm(i),...
          trx(fly1).a_mm(i),trx(fly1).b_mm(i),trx(fly1).theta_mm(i),...
          xnose2,ynose2);
      end
      
      % angle of fly1's vision subtended by fly2
      for t = t0:t1,
        i = t - offi;
        j = t - offj;
        anglesub(fly2,i) = anglesubtended(...
          trx(fly1).x_mm(i),trx(fly1).y_mm(i),trx(fly1).a_mm(i),trx(fly1).b_mm(i),trx(fly1).theta_mm(i),...
          trx(fly2).x_mm(j),trx(fly2).y_mm(j),trx(fly2).a_mm(j),trx(fly2).b_mm(j),trx(fly2).theta_mm(j),...
          params.fov);
      end
      
      % velocity difference
      magveldiff(fly2,i0:i1-1) = ...
        sqrt( (diff(trx(fly1).x_mm(i0:i1))-diff(trx(fly2).x_mm(j0:j1))).^2 + ...
        (diff(trx(fly1).y_mm(i0:i1))-diff(trx(fly2).y_mm(j0:j1))).^2 )./trx(fly1).dt;
      if i1 < trx(fly1).nframes && i1 -1 >= i0,
        magveldiff(fly2,i1) = magveldiff(fly2,i1-1);
      end
      
      % velocity in direction of other fly
      if i1 < trx(fly1).nframes,
        if i1 - i0 + 1 > 0,
          veltoward(fly2,i0:i1) = (dx(1:end).*diff(trx(fly1).x_mm(i0:i1+1)) + ...
            dy(1:end).*diff(trx(fly1).y_mm(i0:i1+1)))./trx(fly1).dt;
        end
      else
        if i1 - i0 > 0,
          veltoward(fly2,i0:i1-1) = (dx(1:end-1).*diff(trx(fly1).x_mm(i0:i1)) + ...
            dy(1:end-1).*diff(trx(fly1).y_mm(i0:i1)))./trx(fly1).dt;
        end
      end
      
      % orientation of fly2 relative to orientation of fly1
      absthetadiff(fly2,i0:i1) = abs(modrange(trx(fly2).theta_mm(j0:j1) - trx(fly1).theta_mm(i0:i1),-pi,pi));
      
      % velocity direction of fly2 relative to fly1's velocity direction
      absphidiff(fly2,i0:i1-1) = abs(modrange(trx(fly2).phi(j0:j1-1)-trx(fly1).phi(i0:i1-1),-pi,pi));
      if i1 < trx(fly1).nframes && i1 -1 >= i0,
        absphidiff(fly2,i1) = absphidiff(fly2,i1-1);
      end
      
      % direction to fly2 from fly1
      absanglefrom1to2(fly2,i0:i1) = abs(modrange(atan2(dy,dx)-trx(fly1).theta_mm(i0:i1),-pi,pi));
      
    end % end loop over fly2
    
    % closest fly according to centroid distance
    [closestfly(fly1).dcenter,closestfly(fly1).closestfly_center] = min(dcenter,[],1);
    units.dcenter = parseunits('mm');
    units.closestfly_center = parseunits('unit');
    
    % closest fly according to dnose2ell
    [closestfly(fly1).dnose2ell,closestfly(fly1).closestfly_nose2ell] = min(dnose2ell,[],1);
    units.dnose2ell = parseunits('mm');
    units.closestfly_nose2ell = parseunits('unit');
    
    % closest fly according to dell2nose
    [closestfly(fly1).dell2nose,closestfly(fly1).closestfly_ell2nose] = min(dell2nose,[],1);
    units.dell2nose = parseunits('mm');
    units.closestfly_ell2nose = parseunits('unit');
    
    % closest fly according to angle subtended
    [closestfly(fly1).anglesub,closestfly(fly1).closestfly_anglesub] = max(anglesub,[],1);
    units.anglesub = parseunits('rad');
    units.closestfly_anglesub = parseunits('unit');
    
    % magveldiff for each of these measured closest flies
    idx = 1:trx(fly1).nframes-1;
    closestfly(fly1).magveldiff_center = magveldiff(sub2ind(size(magveldiff),closestfly(fly1).closestfly_center(idx),idx));
    units.magveldiff_center = parseunits('mm/s');
    closestfly(fly1).magveldiff_nose2ell = magveldiff(sub2ind(size(magveldiff),closestfly(fly1).closestfly_nose2ell(idx),idx));
    units.magveldiff_nose2ell = parseunits('mm/s');
    closestfly(fly1).magveldiff_ell2nose = magveldiff(sub2ind(size(magveldiff),closestfly(fly1).closestfly_ell2nose(idx),idx));
    units.magveldiff_ell2nose = parseunits('mm/s');
    closestfly(fly1).magveldiff_anglesub = magveldiff(sub2ind(size(magveldiff),closestfly(fly1).closestfly_anglesub(idx),idx));
    units.magveldiff_anglesub = parseunits('mm/s');
    
    % veltoward for each of these measured closest flies
    idx = 1:trx(fly1).nframes-1;
    closestfly(fly1).veltoward_center = veltoward(sub2ind(size(veltoward),closestfly(fly1).closestfly_center(idx),idx));
    units.veltoward_center = parseunits('mm/s');
    closestfly(fly1).veltoward_nose2ell = veltoward(sub2ind(size(veltoward),closestfly(fly1).closestfly_nose2ell(idx),idx));
    units.veltoward_nose2ell = parseunits('mm/s');
    closestfly(fly1).veltoward_ell2nose = veltoward(sub2ind(size(veltoward),closestfly(fly1).closestfly_ell2nose(idx),idx));
    units.veltoward_ell2nose = parseunits('mm/s');
    closestfly(fly1).veltoward_anglesub = veltoward(sub2ind(size(veltoward),closestfly(fly1).closestfly_anglesub(idx),idx));
    units.veltoward_anglesub = parseunits('mm/s');
    
    % absthetadiff for each of these measured closest flies
    idx = 1:trx(fly1).nframes;
    closestfly(fly1).absthetadiff_center = absthetadiff(sub2ind(size(absthetadiff),closestfly(fly1).closestfly_center(idx),idx));
    units.absthetadiff_center = parseunits('rad');
    closestfly(fly1).absthetadiff_nose2ell = absthetadiff(sub2ind(size(absthetadiff),closestfly(fly1).closestfly_nose2ell(idx),idx));
    units.absthetadiff_nose2ell = parseunits('rad');
    closestfly(fly1).absthetadiff_ell2nose = absthetadiff(sub2ind(size(absthetadiff),closestfly(fly1).closestfly_ell2nose(idx),idx));
    units.absthetadiff_ell2nose = parseunits('rad');
    closestfly(fly1).absthetadiff_anglesub = absthetadiff(sub2ind(size(absthetadiff),closestfly(fly1).closestfly_anglesub(idx),idx));
    units.absthetadiff_anglesub = parseunits('rad');
    
    % absphidiff for each of these measured closest flies
    idx = 1:trx(fly1).nframes-1;
    closestfly(fly1).absphidiff_center = absphidiff(sub2ind(size(absphidiff),closestfly(fly1).closestfly_center(idx),idx));
    units.absphidiff_center = parseunits('rad');
    closestfly(fly1).absphidiff_nose2ell = absphidiff(sub2ind(size(absphidiff),closestfly(fly1).closestfly_nose2ell(idx),idx));
    units.absphidiff_nose2ell = parseunits('rad');
    closestfly(fly1).absphidiff_ell2nose = absphidiff(sub2ind(size(absphidiff),closestfly(fly1).closestfly_ell2nose(idx),idx));
    units.absphidiff_ell2nose = parseunits('rad');
    closestfly(fly1).absphidiff_anglesub = absphidiff(sub2ind(size(absphidiff),closestfly(fly1).closestfly_anglesub(idx),idx));
    units.absphidiff_anglesub = parseunits('rad');
    
    % anglefrom1to2 for each of these measured closest flies
    idx = 1:trx(fly1).nframes;
    closestfly(fly1).absanglefrom1to2_center = absanglefrom1to2(sub2ind(size(absanglefrom1to2),closestfly(fly1).closestfly_center(idx),idx));
    units.absanglefrom1to2_center = parseunits('rad');
    closestfly(fly1).absanglefrom1to2_nose2ell = absanglefrom1to2(sub2ind(size(absanglefrom1to2),closestfly(fly1).closestfly_nose2ell(idx),idx));
    units.absanglefrom1to2_nose2ell = parseunits('rad');
    closestfly(fly1).absanglefrom1to2_ell2nose = absanglefrom1to2(sub2ind(size(absanglefrom1to2),closestfly(fly1).closestfly_ell2nose(idx),idx));
    units.absanglefrom1to2_ell2nose = parseunits('rad');
    closestfly(fly1).absanglefrom1to2_anglesub = absanglefrom1to2(sub2ind(size(absanglefrom1to2),closestfly(fly1).closestfly_anglesub(idx),idx));
    units.absanglefrom1to2_anglesub = parseunits('rad');
    
    % change in various parameters
    units.ddcenter = parseunits('mm/s');
    closestfly(fly1).ddnose2ell = diff(closestfly(fly1).dnose2ell)./trx(fly1).dt;
    units.ddnose2ell = parseunits('mm/s');
    closestfly(fly1).ddell2nose = diff(closestfly(fly1).dell2nose)./trx(fly1).dt;
    units.ddell2nose = parseunits('mm/s');
    closestfly(fly1).danglesub = diff(closestfly(fly1).anglesub)./trx(fly1).dt;
    units.danglesub = parseunits('rad/s');
    
  end
  
  if isfile,
    fprintf('Saving closest fly measurements to file %s\n',closestfly_file);
    save(closestfly_file,'closestfly','units');
  end
end

% add these fields to trx
for i = 1:length(fns),
  fn = fns{i};
  for fly = 1:nflies,
    trx(fly).(fn) = closestfly(fly).(fn);
  end
end
