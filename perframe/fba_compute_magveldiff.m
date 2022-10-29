% magnitude of difference in velocity between fly and closest fly based on
% type. 
function [data,units] = compute_magveldiff(trx,n,type)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
OLDVERSION = false;

for i1 = 1:nflies,

  fly1 = flies(i1);
  
  % allocate
  data{i1} = nan(1,trx(fly1).nframes-1);
  
  % fly closest to fly1 according to type
  closestfly = trx(fly1).(['closestfly_',type]);
  
  % velocity of fly1
  dx1 = diff(trx(fly1).x_mm,1,2);
  dy1 = diff(trx(fly1).y_mm,1,2);

  % loop over all flies
  for i2 = 1:nflies,
    
    fly2 = flies(i2);
    if i1 == i2, continue; end
    
    % frames where this fly is closest
    idx = find(closestfly(1:end-1) == fly2);
    
    if isempty(idx), continue; end
    
    off = trx(fly1).firstframe - trx(fly2).firstframe;
    % be careful using the last frame of fly2
    usinglastframe = idx(end)+off==trx(fly2).nframes;      
      
    % velocity from these frames to the next frame
    if usinglastframe,
      
      % if this has the last frame, use the previous frame's velocity
      dx2 = trx(fly2).x_mm(off+idx(1:end-1)+1)-trx(fly2).x_mm(off+idx(1:end-1));
      dy2 = trx(fly2).y_mm(off+idx(1:end-1)+1)-trx(fly2).y_mm(off+idx(1:end-1));
      if trx(fly2).nframes == 1,
        dx2(end+1) = 0; %#ok<AGROW>
        dy2(end+1) = 0; %#ok<AGROW>
      else
        dx2(end+1) = trx(fly2).x_mm(end)-trx(fly2).x_mm(end-1); %#ok<AGROW>
        dy2(end+1) = trx(fly2).y_mm(end)-trx(fly2).y_mm(end-1); %#ok<AGROW>
      end
    else
      dx2 = trx(fly2).x_mm(off+idx+1)-trx(fly2).x_mm(off+idx);
      dy2 = trx(fly2).y_mm(off+idx+1)-trx(fly2).y_mm(off+idx);
    end

    % magnitude of difference between velocity vectors
    data{i1}(idx) = sqrt((dx1(idx)-dx2).^2 + (dy1(idx)-dy2).^2);
  end
  if ~OLDVERSION,
    data{i1} = data{i1}./trx(fly1).dt;
  end
end

units = parseunits('mm/s');