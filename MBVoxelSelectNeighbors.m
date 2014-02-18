function [D,idxremove] = MBVoxelSelectNeighbors(pos,maxdneighbor,areaopen,allowdistantneighbors)

maxpos = max(round(pos),[],1);

isdata = false(maxpos);
isdata(sub2ind(maxpos,round(pos(:,1)),round(pos(:,2)),round(pos(:,3)))) = 1;

% remove small connected components
isdata_open = bwareaopen(isdata,areaopen);
isdata_removed = isdata & ~isdata_open;
[xremove,yremove,zremove] = ind2sub(maxpos,find(isdata_removed));

idxremove = find( ismember(round(pos),[xremove,yremove,zremove],'rows') );

pos(idxremove,:) = [];

ndata = size(pos,1);

% to find the neighbors of a given voxel u, we approximately solve the set
% cover problem, where one voxel v covers another voxel w if the distance from
% w to the line (u,v) is <= sqrt(.5^2*3)

% speedups:
% (exact) only consider voxels on the perimeters of the connected
% components
% (approximate) go through the possible neighbors in order of distance to
% the voxel, and greedily add

% also, neighbor relationships must be symmetric, so we add neighbors in
% both directions

% pixel distance in locations to neighbors
D = inf(ndata,ndata,'single');

% set to 1 for direct neighbors (voxels such that if we round their
% positions they will be  of distance <= 1)
for i = 1:ndata,
  pos1 = pos(i,:);
  
  for dir = 1:3,
    
    idxnear = true(ndata,1);
    idxnear(i) = false;
    for j = 1:3,
      if j == dir,
        idxnear = idxnear & abs(pos(:,j)-pos1(j)) < 1.5;
      else
        idxnear = idxnear & abs(pos(:,j)-pos1(j)) < .5;
      end
    end
    dcurr = sum(abs(bsxfun(@minus,pos(idxnear,:),pos1)),2);
    D(i,idxnear) = dcurr;
    D(idxnear,i) = dcurr;
    
  end
end


% % find voxels that do not have all their neighbors (on the perimeter)
% isperim = bwperim(isdata,6);
% 
% % x,y,z locations of these perimeter pixels
% [xperim,yperim,zperim] = ind2sub(maxpos,find(isperim));
% % data index for these perimeter pixels
% idxperim = dataidx(isperim);

% find neighbors to add

% maximum allowed distance from line to a neighbor to consider that line
% occluded

clf;

plot3(pos(:,1),pos(:,2),pos(:,3),'k.');
hold on;
for i = 1:ndata,
  idx = ~isinf(D(i,:));
  z = zeros(nnz(idx),1);
  plot3([pos(i,1)+z,pos(idx,1)]',[pos(i,2)+z,pos(idx,2)]',[pos(i,3)+z,pos(idx,3)]','-');
end
title('initial neighbors');
drawnow;

if allowdistantneighbors,
  maxd = sqrt(.5^2*3);
  
  % loop through all the voxels that do not have all their neighbors yet
  for i = 1:ndata,
    
    % voxel at pos(i,1),pos(i,2),pos(i,3)
    x = pos(i,1); y = pos(i,2); z = pos(i,3);
    poscurr = pos(i,:);
    
    % grab a box around this voxel
    off1 = poscurr-maxdneighbor;
    off2 = poscurr+maxdneighbor;
    
    % voxels in this box
    idxnear = all(bsxfun(@ge,pos,off1) & bsxfun(@le,pos,off2),2);
    idxnear(i) = false;
    idxnear = find(idxnear);

    % go through all voxels nearby in order of their distance to the current voxel
    distmaybenear = sum(bsxfun(@minus,pos(idxnear,:),poscurr).^2,2);
    [~,order] = sort(distmaybenear);
    xmaybenear = pos(idxnear,1);
    ymaybenear = pos(idxnear,2);
    zmaybenear = pos(idxnear,3);
    
    % initialize neighbors to include neighbors that have already been
    % added
    tmp = find(~isinf(D(i,:)));
    xdatanear = pos(tmp,1);
    ydatanear = pos(tmp,2);
    zdatanear = pos(tmp,3);
    
    fprintf('%d/%d: %d neighbors initially ... ',i,numel(ndata),numel(xdatanear));
    
    for j = order(:)',
%       if distmaybenear(j) <= 1,
%         % these should already be in the initial datanear
%         if ~any(xdatanear==xmaybenear(j)&ydatanear==ymaybenear(j)&zdatanear==zmaybenear(j)),
%           error('Sanity check: adjacent voxel not yet in initial neighbor set');
%         end
%         continue;
%       end
      
      % already in the neighbor set
      if any(xdatanear==xmaybenear(j)&ydatanear==ymaybenear(j)&zdatanear==zmaybenear(j)),
        continue;
      end
      
      % if no voxels near yet
      if isempty(xdatanear),
        dlineseg = [];
      else
        % compute the distance from all neighbors to the line joining these
        % voxels
        p = [xdatanear';ydatanear';zdatanear'];
        a = repmat(poscurr',[1,size(p,2)]);
        b = repmat([xmaybenear(j);ymaybenear(j);zmaybenear(j)],[1,size(p,2)]);
        dlineseg = distancePointToLineSegmentV(p,a,b);
      end
      
      % if no neighbors are occluding this line yet add it to the neighbor
      % set
      if ~any(dlineseg<=maxd),
        xdatanear(end+1,1) = xmaybenear(j); %#ok<AGROW>
        ydatanear(end+1,1) = ymaybenear(j); %#ok<AGROW>
        zdatanear(end+1,1) = zmaybenear(j); %#ok<AGROW>
        % index into data
        jj = idxnear(j);
        if ~isinf(D(i,jj)),
          error('Sanity check: distance already set for these voxels');
        end
        % add it to distance matrix
        D(i,jj) = abs(xmaybenear(j)-x)+...
          abs(ymaybenear(j)-y)+...
          abs(zmaybenear(j)-z);
        D(jj,i) = D(i,jj);
      end
      
    end
    
    fprintf('total %d neighbors after adding.\n',numel(xdatanear));
    
    
  end
  
end
